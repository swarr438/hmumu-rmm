'''Transform RMMs(7x7) into RGB images(3x5x5).'''

from math import log10, floor
from pathlib import Path
from PIL import Image
import numpy as np
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)
# print(os.getcwd())

# rgb algorithm
def find_exp(number) -> int:
    base10 = log10(abs(number))
    return abs(floor(base10))

def log2r01(number):
    base10 = 1
    if number>0:
        base10 = abs(log10(number))
        # if base10<10:
        #     return base10/10
        # else:
        #     return 1
    
    return base10

def num2rgb(number):
    # number must be between 0 and 1
    if number!=0:
        pow = find_exp(number)
        float = log10(abs(number))-pow
        gb = floor(float*10**5)# (pow+4)
        if gb>500000:
            r = pow*30-15
        else:
            r= pow*30
        if r>255:
            r=255
        g = gb//256
        b = gb%256
        return r,g,b
    else:
        return 255,255,255

def num2rgb2(number):
    if number!=0:
        pow = find_exp(number)
        float = abs(log10(number))
        g = floor((float-pow+1)*255)
        r = (10-pow)*25
        if r<0:
            r=0
        return [r,g,0]
    else:
        return [0,0,0]

filenames = ["dataset0"] #["dataset0","dataset1","dataset2"]

for filename in filenames:
    print("-------------------------")
    print(f"Generating {filename}...")
    Path(f"{dir_path}\{filename}").mkdir(parents=False, exist_ok=True)
    count = 0
    index = open(filename+".csv","w")
    index.write('input,output\n')
    with open(filename+".txt","r") as ifile:
        for line in ifile:
            count+=1
            ls = line.split(" ")
            ncells = len(ls)-2
            nin = floor(np.sqrt(ncells))
            # arr = np.zeros((nin,nin,3))
            # arrgr = np.zeros((nin,nin))
            
            # Rollback to T2N2
            shapedls = np.asarray(ls[1:ncells+1],dtype=float).reshape(nin,nin)
            if float(ls[1])<0.5:
                arr = np.delete(shapedls,[3,4],0)
                arr = np.delete(arr,[3,4],1)
                b = np.zeros(arr.shape)
            elif float(ls[1])>0.5:
                arr = np.delete(shapedls,[1,2],0)
                arr = np.delete(arr,[1,2],1)
                b = np.zeros(arr.shape)
                b[:,1]=np.ones(len(b[:,1]))*255
                b[:,2]=np.ones(len(b[:,2]))*255
                b[1,:]=np.ones(len(b[:,1]))*255
                b[2,:]=np.ones(len(b[:,2]))*255
            else:
                arr = np.delete(shapedls,[2,4],0)
                arr = np.delete(arr,[2,4],1)
                b = np.zeros(arr.shape)
                b[:,2]=np.ones(len(b[:,2]))*255
                b[2,:]=np.ones(len(b[:,2]))*255
            
            # print(arr)
            # imgg = Image.fromarray(arr*255).convert('L')
            # imgg = imgg.resize((500, 500),resample=PIL.Image.Resampling.NEAREST)
            # imgg.show()
            # imgg.save('gray.png')

            # delete diagonal values
            np.fill_diagonal(arr,0)
            np.fill_diagonal(b,0)
            arrgb = np.zeros((5,5,3))
            for i in range(arr.shape[0]):
                for j in range(arr.shape[1]):
                    arrgb[i][j] = num2rgb2(arr[i][j])
            
            arrgb[:,:,2] = b

            # for i in range(nin):
            #     for j in range(nin):
            #             pos = 1+i*nin+j
            #             val = float(ls[pos])
            #             arr[i][j][0],arr[i][j][1],arr[i][j][2] = num2rgb(val)
            #             #arrgr[i][j] = log2r01(val)

            img = Image.fromarray(arrgb.astype('uint8'), 'RGB')
            # img = img.resize((500, 500),resample=PIL.Image.Resampling.NEAREST)
            #img.show()
            #img.save('rgb.png')
            imgname = str(count)+".png"
            img.save(f"{filename}\{imgname}")
            index.write(f"{imgname},{ls[0]}\n")

            if count%10000==0:
                print(f"Finished {count} events...")
            
    print(f"Finished {count} events in total")
    index.close()