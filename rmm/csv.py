# generate datasets from channels

channel = {"ttbar":840, "diboson":160, "Z":9000, "ttH":220, "VH":280, "VBF":700, "ggF":8800}
chout = {"ttbar":"0", "diboson":"0", "Z":"0", "ttH":"1", "VH":"2", "VBF":"3", "ggF":"4"}

# num = 20000
fname = "dataset"

for j in range(3):
    ofname = fname+str(j)
    # ofname = 'dataset3'
    ofile = open(f"{ofname}.txt","w")
    # ofile.write("21000 49 2\n")

    for x,y in channel.items():
        with open(x+".txt","r") as ifile:
            lines = ifile.readlines()[j*y:(j+1)*y]
            for i in lines:
                ofile.write(chout[x]+" ")
                ofile.write(i)

    ofile.close()

    with open(f"{ofname}.txt","r") as ofile:
        n = len(ofile.readlines())
        print(f'{ofname}: {n} events...')
