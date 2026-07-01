![jamlogo](logos/jam.jpg)

# THIS LIBRARY IS NOW DEPRECATED - PLEASE VISIT https://colab.research.google.com/github/pitonyak25/jam3d_dev_lib/blob/main/JAM3D_Library.ipynb

# JAM3DLIB

## List of Fits

|Reference        |   tag    |  Data Included                                               |   Functions Extracted       |
|-----------------|----------|--------------------------------------------------------------|-----------------------------|
|arXiv:2002.08384 |JAM3D_2020|SIDIS (Sivers, Collins), SIA (Collins), DY (Sivers), pp (AN)  |Sivers, transversity, Collins|                                                                      

## Installation

### clone the repos

 ```git clone  git@github.com:QCDHUB/jam3dlib.git```

 ```git clone  git@github.com:QCDHUB/jam3d.git```

### prepare a python 2.7 environment 

- Download and install anaconda3 (python 3.7) into your system from https://www.anaconda.com/distribution/ (python 2.7 is no longer supported but is needed to run jam3d)

- Create and activate a python 2.7 enviroment

  ```conda create -n python2 python=2.7 anaconda```
  
  ```conda activate python2```

- Install subprocess32

  ```cd jam3dlib```

  ```conda install subprocess32```

- Install all the dependencies

  ```pip install -r dependencies```


## Getting started

- Change into the jam3d directory 

  ```cd ../jam3d```

- Source the setup

  ```source setup.bash```

- Change back into the jam3dlib directory

  ```cd ../jam3dlib```

- We have provided two files in the jam3dlib repo with simple
examples to evaluate the TMD and CT3 PDFs/FFs and asymmetries

   - example.ipynb: detailed example jupyter notebook

     ```jupyter notebook example.ipynb```

   - example.py: bare-bones template python file

     ```python example.py```



## Authors (Please contact Daniel, Alexei, or Nobuo with questions)

- Justin Cammarota

- Leonard Gamberg

- Zhong-Bo Kang

- Joshua A. Miller

- Daniel Pitonyak (pitonyak@lvc.edu)

- Alexei Prokudin (prokudin@jlab.org)

- Ted C. Rogers

- Nobuo Sato (nsato@jlab.org)



## Institutions

![l](logos/LVC.jpeg)
![logo](logos/PSU_BKO_RGB_2C.png)
![logo](logos/odu.png)
![logo](logos/JLab_logo_white1.jpg)
![logo](logos/NSF_4-Color_bitmap_Logo.png)
![logo](logos/RGB_Color-Seal_Green-Mark_SC_Horizontal.png)
