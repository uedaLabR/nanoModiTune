import filter.NNModel as NNModel

model = NNModel.getModel()
model.summary()
model.save_weights('/mnt/share/ueda/RNA004/resource/test.weights.h5')

newmodel = NNModel.getModel()
newmodel.load_weights('/mnt/share/ueda/RNA004/resource/test.weights.h5')
newmodel.summary()

