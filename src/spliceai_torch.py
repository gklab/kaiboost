"""Run the official SpliceAI Keras .h5 weights with PyTorch (no TensorFlow).

TensorFlow crashes at import in this conda env (libc++ mutex abort), so the
five bundled models are executed by walking the Keras Functional graph stored
in each .h5 `model_config` attribute. Layer types present: InputLayer, Conv1D
(linear/softmax), BatchNormalization, Activation(relu), Add, Cropping1D.

The .predict(x) contract matches what spliceai.utils.get_delta_scores expects:
numpy (B, L, 4) channels-last in -> numpy (B, L-10000, 3) softmax out.
"""
import json

import h5py
import numpy as np
import torch
import torch.nn.functional as F


class TorchSpliceAI:
    def __init__(self, h5path):
        with h5py.File(h5path, "r") as f:
            cfg = json.loads(f.attrs["model_config"])["config"]
            self.layers = cfg["layers"]
            self.output = cfg["output_layers"][0][0]
            self.params = {}
            mw = f["model_weights"]
            for lname in mw:
                grp = mw[lname]
                names = [n.decode() if isinstance(n, bytes) else n
                         for n in grp.attrs.get("weight_names", [])]
                store = {}
                for wn in names:
                    store[wn.split("/")[-1]] = np.array(grp[wn])
                if store:
                    self.params[lname] = store

        self.tensors = {}
        for lname, ws in self.params.items():
            t = {}
            for k, v in ws.items():
                arr = v.astype(np.float32)
                if k == "kernel:0":            # keras (W, Cin, Cout) -> torch (Cout, Cin, W)
                    arr = np.ascontiguousarray(arr.transpose(2, 1, 0))
                t[k] = torch.from_numpy(arr)
            self.tensors[lname] = t

    def _conv(self, x, layer):
        c = layer["config"]
        w = self.tensors[layer["name"]]["kernel:0"]
        b = self.tensors[layer["name"]].get("bias:0")
        dil = c["dilation_rate"][0]
        ksz = c["kernel_size"][0]
        if c["padding"] == "same" and ksz > 1:
            total = dil * (ksz - 1)
            x = F.pad(x, (total // 2, total - total // 2))
        y = F.conv1d(x, w, b, dilation=dil)
        if c["activation"] == "softmax":
            y = F.softmax(y, dim=1)
        elif c["activation"] == "relu":
            y = F.relu(y)
        return y

    def _bn(self, x, layer):
        ws = self.tensors[layer["name"]]
        eps = layer["config"].get("epsilon", 1e-3)
        gamma, beta = ws["gamma:0"], ws["beta:0"]
        mean, var = ws["moving_mean:0"], ws["moving_variance:0"]
        shape = (1, -1, 1)
        return (x - mean.view(shape)) / torch.sqrt(var.view(shape) + eps) \
            * gamma.view(shape) + beta.view(shape)

    @torch.no_grad()
    def predict(self, x, batch_size=32):
        x = np.ascontiguousarray(x, dtype=np.float32)
        t = torch.from_numpy(x).permute(0, 2, 1)      # (B, 4, L)
        acts = {}
        for layer in self.layers:                     # keras config is topological
            cls, name = layer["class_name"], layer["name"]
            if cls == "InputLayer":
                acts[name] = t
                continue
            srcs = [acts[n[0]] for n in layer["inbound_nodes"][0]]
            if cls == "Conv1D":
                acts[name] = self._conv(srcs[0], layer)
            elif cls == "BatchNormalization":
                acts[name] = self._bn(srcs[0], layer)
            elif cls == "Activation":
                acts[name] = F.relu(srcs[0])
            elif cls == "Add":
                acts[name] = sum(srcs)
            elif cls == "Cropping1D":
                a, b = layer["config"]["cropping"]
                acts[name] = srcs[0][:, :, a:srcs[0].shape[2] - b]
            else:
                raise ValueError(f"unhandled layer {cls}")
        return acts[self.output].permute(0, 2, 1).numpy()
