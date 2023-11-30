# hmumu-rmm

This project, inspired by S. Chekanov's [Map2RMM](https://github.com/chekanov/Map2RMM), aims to search for the SM Higgs to dimuon decay by combining rapidity-mass matrices and convolutional neural network. 

1. To begin, run **rmm/rmmf.C** in the shell with
```
root -q rmmf.C
```
and input the channel name. There are seven channels available: _ttbar_, _diboson_, _Z_, _ttH_, _VH_, _VBF_, _ggF_. The corresponding Monte Carlo file ids and event ids can be adjusted in the **{channel}/index/** directory. In this repository we only give an example of the _ggF_ channel.
In order to complete the training dataset, a loop over all channels is required and finally you will get a _txt_ file that store the RMM values of events for each channel.

>[!NOTE]
>The input data are stored on the CERN server so you must run it, for example, on LXPLUS.

2. Next, run **rmm/csv.py**. This will select events from each channel and generate the datasets for training.

3. Third, run **rmm2rgb.py**, which visualizes all events in the datasets. The output images will be used as the input to CNN.

4. Now you can start training. Since GPU is necessary, we recommend using the colab notebook.
