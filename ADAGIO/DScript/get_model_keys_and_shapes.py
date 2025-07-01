from safetensors.torch import load_file

path = "model.safetensors"
state_dict = load_file(path)

for key, tensor in state_dict.items():
    print(f"{key}: {tuple(tensor.shape)}")
