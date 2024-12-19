import filter.NNModel as NNModel

# model = NNModel.getModel()
# model.summary()
# model.save_weights('/mnt/share/ueda/RNA004/resource/test_stats.weights.h5')
#
# newmodel = NNModel.getModel()
# newmodel.load_weights('/mnt/share/ueda/RNA004/resource/test_stats.weights.h5')
# newmodel.summary()


import numpy as np

def average_by_group(array, group_size=3):

    array = np.array(array)  # Ensure input is a NumPy array
    truncated_length = (len(array) // group_size) * group_size  # Compute truncation length
    truncated_array = array[:truncated_length]  # Truncate to a multiple of group_size
    reshaped = truncated_array.reshape(-1, group_size)  # Reshape into groups
    return reshaped.mean(axis=1)  # Compute the mean along the groups

# Example usage
original_array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
averaged_array = average_by_group(original_array, group_size=3)

print("Original array:", original_array)
print("Averaged array:", averaged_array)



# Example usage
original_array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
averaged_array = average_by_group(original_array, group_size=3)

print("Original array:", original_array)
print("Averaged array:", averaged_array)
