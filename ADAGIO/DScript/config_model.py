"""
must run from DIAMONDS25/ADAGIO
"""

import sys
import os

# Add project root to path
sys.path.append(os.path.abspath(os.path.dirname(__file__) + "/.."))

from safetensors.torch import load_file
from DScript.models.interaction import DSCRIPTModel

config = {
    "emb_nin": 6165,
    "emb_nout": 100,
    "emb_dropout": 0.0,
    "con_embed_dim": 200,
    "con_hidden_dim": 50,
    "con_width": 7,
    "use_cuda": False,
    "do_w": True,
    "do_sigmoid": True,
    "do_pool": False,
    "pool_size": 9,
    "theta_init": 1,
    "lambda_init": 0,
    "gamma_init": 0,
}

# Reconstruct model
model = DSCRIPTModel(**config)
print("Model conv input channels:", model.contact.hidden.conv.weight.shape[1])


# Load weights
state_dict = load_file("DScript/model.safetensors")
model.load_state_dict(state_dict)

# Save properly
model.save_pretrained("DScript/saved_model_dir")