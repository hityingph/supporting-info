from dpdata import LabeledSystem,MultiSystems


train_multi_systems = MultiSystems.from_file(
    file_name="train.xyz", fmt="quip/gap/xyz")

train_multi_systems.to_deepmd_raw("./train_dataset")
train_multi_systems.to_deepmd_npy("./train_dataset")


test_multi_systems = MultiSystems.from_file(
    file_name="test.xyz", fmt="quip/gap/xyz")

test_multi_systems.to_deepmd_raw("./test_dataset")
test_multi_systems.to_deepmd_npy("./test_dataset")