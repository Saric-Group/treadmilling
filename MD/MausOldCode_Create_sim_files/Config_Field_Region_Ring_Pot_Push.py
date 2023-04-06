
import numpy as np 
import matplotlib.pyplot as plt
import time
import scipy.spatial as spa
#Crearemos un script que cree el archivo necesario 

Atoms=12000
Fil_L=30
Box=[-100.0,100.0,-150.0,150.0 , -0.3,0.3]
f_swim=50.0
Nsim_ini=1
Nsim_fin=2
FField=10.0
sigma2=25
f_chem=0.0
f_chem2=100
def dist(x,y,x2,y2,Box):
    Lx=Box[1]-Box[0]
    Ly=Box[3]-Box[2]
    
    #IN this case you can have PBC in X and y
    dx=min( abs(x-Lx-x2),abs(x+Lx-x2),abs(x-x2))
    dy=min( abs(y-Ly-y2),abs(y+Ly-y2),abs( y-y2))

    return np.sqrt(dx*dx +dy*dy ) 



def dist_nopbc(x,y,x2,y2):
    return np.sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2))

def Atom_pos_Uniform(Atoms, Fil_L, Box):

    Mols=int(Atoms/Fil_L)
    X=[]
    Y=[]
    Z=[]
    Lx=Box[1]-Box[0]
    Ly=Box[3]-Box[2]
    print('The length of the box is {}X{}'.format(Lx,Ly))
    #We have a number of atoms that would like to be distributed uniformly.
    #We could calculate the distribution of heads
    #Since the field aligns them with the x axis lets have them pointing in the y axis
    #We need to find how many rows we can put
    #W
    extra_height=12
    height=(Fil_L+extra_height)
    Nrows=int(Ly/height)
    print('The number of rows is {}'.format(Nrows))
    #This will tell me how much rows y can make
    #to find the number of colums i have to do 
    Ncolumns= int(Mols/Nrows)+1
    print('The number of columns is {}'.format(Ncolumns))
    
    width=Lx/Ncolumns
    print('The width of each section is {}, and its height is {}'.format(width,height))

    Actual_row=0
    Actual_column=0
    Origin_x=Box[0]
    Origin_y=Box[2]
    for i in range(Mols):
        sign=1
        x_head= width/2+width*Actual_column+Origin_x
        if(i%2==0):
            #MOlecule pointing upwards
            y_head= Actual_row*height + extra_height/2+Origin_y
            sign=1
        else:
            #Molecule pointing downwards
            y_head= Actual_row*height + extra_height/2+Fil_L+Origin_y
            sign=-1
        #Now i have the heads I have to add the body

        y_mol=y_head +sign*np.array(range(Fil_L))
        y_mol=y_mol.tolist()
        x_mol=np.ones_like(y_mol)*x_head 
        x_mol= x_mol.tolist()
        #NOw i have the x_mol and y_mol for each one
        X=X+x_mol
        Y=Y+y_mol
        plt.plot(x_mol,y_mol)

        #Now we need to advance in column and row
        Actual_column+=1
        if(Actual_column==Ncolumns):
            Actual_row+=1
            Actual_column=0
        


    plt.ylim(Box[2],Box[3])
    plt.xlim(Box[0],Box[1])
    # plt.show()


    Z=[0]*len(X)

    return [X,Y,Z]

# start=time.time()
# Atom_pos_Uniform(Atoms,Fil_L,Box)
# end=time.time()
# print('This operation for {} atoms took {} seconds'.format(Atoms,end-start))

def Atom_pos(Atoms, Fil_L, Box,g):
    rejects=0
    #THe idea is to create particles within the box
    Mols=int(Atoms/Fil_L)
    X=[]
    Y=[]
    Z=[]
    Lx=Box[1]-Box[0]
    Ly=Box[3]-Box[2]
    for i in range(Mols):
        #THis will be for each molecule
        warning=False
        accept=False
        while(not(accept) ):
            if(warning):
                # print('Y se volvio hasta aca')
                warning=False
            x_mol=[]
            y_mol=[]
            x_cand=np.random.random()*(Lx) -Lx/2
            y_cand=np.random.random()*(Ly) -Ly/2
            
            #I will try to find it until it works haha
            for j in range(len(X)):
                distance=dist(x_cand,y_cand,X[j],Y[j],Box)
                if(distance<1.1):
                    continue
            #If you got here you are not close to any rod
            #Now we pick an orientation

            angle=np.random.random()*2*np.pi
            angle=np.random.random()*2*np.pi/6 -np.pi/12 +np.pi/2
            coin=np.random.random()
            if(coin>=0.5):
                angle=angle+np.pi
            # print('{} EL angulo cambia ves \n'.format(angle))
            
            for k in range(Fil_L):
                #Wee start adding the particles
                x_next= x_cand+0.9*np.cos(angle)*k
                y_next= y_cand+0.9*np.sin(angle)*k
                if(x_next<Box[0]):
                    x_next=x_next+Lx
                if(x_next>Box[1]):
                    x_next=x_next-Lx
                if(y_next<Box[2]):
                    y_next=y_next+Ly
                if(y_next>Box[3]):
                    y_next=y_next-Ly
                
                for m in range(len(X)):
                    distance=dist(x_next,y_next,X[m],Y[m],Box)
                    if(distance<1.1):
                        warning=True
                        rejects+=1
                        # print('Se rechaza en i= {}\n'.format(m))
                        continue
            
                x_mol.append(x_next)
                y_mol.append(y_next)
            if(warning):
                rejects+=1
                continue
            # plt.scatter(x_mol,y_mol)
            
                #do we check everytime?
            accept=True
        X=X+x_mol 
        Y=Y+y_mol
        Z=[0]*len(X)
    
    # plt.show()
    # plt.scatter(X,Y)
    # plt.show()
    String='The number of rejects is {}\n'.format(rejects)
    print(String)
    g.write(String)
    return [X,Y,Z]






def Atom_pos_opt(Atoms, Fil_L, Box,g):
    rejects=0
    #THe idea is to create particles within the box
    Mols=int(Atoms/Fil_L)
    X=[]
    Y=[]
    Z=[]
    Lx=Box[1]-Box[0]
    Ly=Box[3]-Box[2]
    print('THe occupancy is {}'.format(Atoms/(Lx*Ly)))
    for i in range(Mols):
        #THis will be for each molecule
        warning=False
        accept=False
    
        Tree=spa.cKDTree(np.array([X,Y]).transpose(),leafsize=16, compact_nodes=True, copy_data=False, balanced_tree=True, boxsize=[Lx,Ly])


        while(not(accept) ):
            if(warning):
                # print('Y se volvio hasta aca')
                warning=False
            x_mol=[]
            y_mol=[]
            x_cand=np.random.random()*(Lx) 
            y_cand=np.random.random()*(Ly) 
            #
            #I have to find the neighbours now
            distance=Tree.query(np.array([x_cand,y_cand] ))
            
            if(distance[0]<1.1):
                warning=True
                continue
           
            #If you got here you are not close to any rod
            #Now we pick an orientation

            angle=np.random.random()*2*np.pi
            # angle=np.random.random()*2*np.pi/6 -np.pi/12 +np.pi/
            # coin=np.random.random()
            # if(coin>=0.5):
            #     angle=angle+np.pi
            # print('{} EL angulo cambia ves \n'.format(angle))
            
            for k in range(Fil_L):
                #Wee start adding the particles
                x_next= x_cand+0.9*np.cos(angle)*k
                y_next= y_cand+0.9*np.sin(angle)*k
                if(x_next<0):
                    x_next=x_next+Lx
                if(x_next>Lx):
                    x_next=x_next-Lx
                if(y_next<0):
                    y_next=y_next+Ly
                if(y_next>Ly):
                    y_next=y_next-Ly
                
                #I need to make sure the new point doesnt crash
                distance=Tree.query(np.array([x_next,y_next] ))
                if(distance[0]<1.1):
                    warning=True
                    continue
           
            
                x_mol.append(x_next)
                y_mol.append(y_next)
            if(warning):
                rejects+=1
                continue
            # plt.scatter(x_mol,y_mol,s=0.5)
            
                #do we check everytime?

            #Here i need to check that the distance between the first and the second is less than 0.9
            distance= dist_nopbc(x_mol[0],y_mol[0],x_mol[1],y_mol[1])
            if(distance>2):
                rejects+=1
                print('Paso este caso de mrd')
                continue
            accept=True

        X=X+x_mol 
        Y=Y+y_mol
    # plt.show()
    X=np.array(X)-Lx/2
    Y=np.array(Y)-Ly/2
    Z=[0]*len(X)
    

    # plt.show()
    # plt.scatter(X,Y)
    # plt.show()
    String='The number of rejects is {}\n'.format(rejects)
    
    print(String)
    g.write(String)
    return [X,Y,Z]





# g=open('./Whatever.txt','w')


# start1=time.time()
# Atom_pos_opt(Atoms,Fil_L,Box,g)
# end1=time.time()
# print('This operation for {} atoms took {} seconds in the optimized way \n'.format(Atoms,end1-start1))



def Input_file(Atoms,Fil_L,Box,Config_dir,f_swim,f_chem,f_chem2,sigma2,FField,Nsim):
    f = open('Inputs/in.MD_RODS_N_{}_L_{}_Lx_{}_Ly_{}_fswim_{}_fchem_{}_fchem2_{}_sigma2_{}_Nsim_{}_FFieldRegion_{}_Pushing'.format(Atoms,Fil_L,Box[1],Box[3],f_swim,f_chem,f_chem2,sigma2,Nsim,FField), 'w')
    Ly=Box[3]-Box[2]

    f.write('units          lj\n')
    f.write('dimension      2\n')
    f.write('\n')
    f.write('atom_style     molecular')
    f.write('\n')
    f.write('read_data      {}'.format(Config_dir))
    f.write('\n')
    f.write('\n')
    f.write('bond_style     harmonic\n')
    f.write('bond_coeff     1 1500 0.9\n')
    f.write('\n')
    f.write('\n')
    f.write('angle_style    harmonic\n')
    f.write('angle_coeff    * 25000 180.0\n')
    f.write('\n')
    f.write('\n')
    f.write('pair_style     lj/cut 1.2\n')
    f.write('pair_coeff     * * 0.01 1 1.1\n')
    f.write('pair_modify    shift true\n')
    f.write('\n')
    f.write('\n')
    f.write('neigh_modify   exclude molecule/intra all\n')
    f.write('\n')
    f.write('\n')
    f.write('timestep 0.001\n')
    f.write('\n')
    f.write('\n')
    f.write('group          head type 1\n')
    f.write('group          neck type 2\n')
    f.write('\n')
    f.write('\n')
    f.write('fix            1 all langevin 1.0 1.0 1.0 134\n')
    f.write('fix            2 all nve\n')
    f.write('fix            3 all enforce2d\n')
    f.write('\n')
    f.write('\n')

    f.write('dump           1 all custom 1000 Dumps/dump_N_{}_L_{}_Lx_{}_Ly_{}_fswim_{}_fchem_{}_fchem2_{}_sigma2_{}_Nsim_{}_FField_{}_pushing.xyz id mol type x y z vx vy vz\n'.format(Atoms,Fil_L,Box[1],Box[3],f_swim,f_chem,f_chem2,sigma2,Nsim,FField))
    f.write('\n')
    f.write('\n')

    f.write('compute        cMol all chunk/atom molecule\n')
   
    f.write('\n')
    f.write('\n')
    # f.write('compute        cCM all com/chunk cMol\n')
    # f.write('fix            Heating all spring/chunk 100 cMol cCM\n')
    # f.write('run            80000\n')
    # f.write('unfix          Heating\n')
    f.write('\n')
    f.write('\n')


    f.write('compute        cCMNeck neck com/chunk cMol\n')
    f.write('compute        cNeck all chunk/spread/atom cMol c_cCMNeck[*]\n')
    f.write('compute        cxu all property/atom xu\n')
    f.write('compute        cyu all property/atom yu\n')
    f.write('compute        cx all property/atom x\n')
    f.write('compute        cy all property/atom y\n')
    f.write('variable       cosalph atom (c_cxu-c_cNeck[1])/((c_cyu-c_cNeck[2])*(c_cyu-c_cNeck[2])+(c_cxu-c_cNeck[1])*(c_cxu-c_cNeck[1]))\n'.format(f_swim))
    f.write('variable       sinalph atom (c_cyu-c_cNeck[2])/((c_cyu-c_cNeck[2])*(c_cyu-c_cNeck[2])+(c_cxu-c_cNeck[1])*(c_cxu-c_cNeck[1]))\n'.format(f_swim))
    f.write('variable       fxHead atom {}*v_cosalph\n'.format(f_swim))
    f.write('variable       fyHead atom {}*v_sinalph\n'.format(f_swim))
    f.write('variable       fxField atom {}*(1.0-v_cosalph)*v_sinalph*v_sinalph*((c_cxu-c_cNeck[1])/abs((c_cxu-c_cNeck[1])))       \n'.format(FField))
    f.write('variable       fyField atom {}*(1.0-v_cosalph)*v_sinalph*-1.0*v_cosalph*((c_cxu-c_cNeck[1])/abs((c_cxu-c_cNeck[1])))  \n'.format(FField))
    f.write('\n')
    f.write('fix            fswimmer head addforce v_fxHead v_fyHead 0\n')
    f.write('\n')
    f.write('\n')
    f.write('region         ForceRegion block {} {} -20.0 20.0 {} {}\n'.format(Box[0],Box[1],Box[4],Box[5]))
    f.write('\n')
    f.write('\n')
    # f.write('fix            fField2 head addforce v_fxField v_fyField 0.0 region ForceRegion\n')
    f.write('fix            fField head addforce v_fxField v_fyField 0.0\n')
    # Tihs is a field that makes particles align with the x axis
    f.write('\n')
    f.write('\n')
    # f.write('variable       fxChem atom {}*v_cosalph*(exp(-1*y*y/2.0)-exp(-1*{}*{}/2.0) )/(1-exp(-1*{}*{}/2.0 ))  \n'.format(f_chem,Ly/2,Ly/2,Ly/2,Ly/2))
    
    f.write('variable       denominator atom 1.0-exp(-1.0*{}/{})\n'.format(Ly*Ly/8,sigma2))
    f.write('variable       numerator   atom exp(-1.0*c_cy*c_cy/(2.0*{}))-exp(-1.0*{}/{})\n '.format(sigma2,Ly*Ly/8,sigma2))
    
    
    f.write('variable       fxChem2 atom {}*(1.0-v_cosalph)*v_sinalph*v_sinalph*((c_cxu-c_cNeck[1])/abs((c_cxu-c_cNeck[1])))*(v_numerator/v_denominator)       \n'.format(f_chem2))
    f.write('variable       fyChem2 atom {}*(1.0-v_cosalph)*v_sinalph*(-1.0)*v_cosalph*((c_cxu-c_cNeck[1])/abs((c_cxu-c_cNeck[1])))*(v_numerator/v_denominator)  \n'.format(f_chem2))
    # THis force is also of alignment but decays exponentially 
    
    f.write('variable       fyChem atom -1*{}*c_cy  \n'.format(f_chem))
    f.write('fix            fChem head addforce 0 v_fyChem 0.0\n')
    # This force is the attractive potential
    
    f.write('fix            fChem2 head addforce v_fxChem2 v_fyChem2 0.0\n')
    
    
    f.write('\n')
    f.write('\n')
    
    f.write('region         PushRegion block -20.0 20.0 -30.0 30.0 {} {}\n'.format(Box[4],Box[5]))
    f.write('run            600000\n')
    f.write(' \n')
    f.write(' \n')
    f.write(' \n')
    f.write(' \n')
    f.write('fix            fPush all addforce 0 5 0 region PushRegion\n')
    f.write('run            100000\n')
    f.write('unfix          fPush')
    f.write(' \n')
    f.write(' \n')

    f.write('run            400000\n')
    # f.close()




    return 


def Config_file(Atoms,Fil_L,Box, Nsim,Bonds, Angles, Angle_types,Atom_types,Bond_types,Config_dir,g):

    # Coords=Atom_pos(Atoms,Fil_L,Box,g)
    Coords=Atom_pos_opt(Atoms,Fil_L,Box,g)
    
    
    #The first part is easy.. you first create de configuration file
    f = open(Config_dir, 'w')
    f.write('First line of this test\n')
    f.write('{} atoms\n'.format(Atoms) )
    f.write('{} bonds\n'.format(Bonds))
    f.write('{} angles\n'.format(Angles))
    f.write('{} atom types\n'.format(Atom_types))
    f.write('{} bond types\n'.format(Bond_types))
    f.write('{} angle types\n'.format(Angle_types))
    f.write('{} {} xlo xhi\n'.format(Box[0],Box[1]))
    f.write('{} {} ylo yhi\n'.format(Box[2],Box[3]))
    f.write('{} {} zlo zhi\n'.format(Box[4],Box[5]))

    #With this the header ends and we start with the more problematic part
    #You need to understand the order 
    #FOr atoms is   id mol type x y z
    f.write('\n')
    f.write('Masses\n')
    f.write('\n')
    
    f.write('1 1\n')
    f.write('2 1\n')
    f.write('3 1\n')

    f.write('\n')
    f.write('Atoms\n')
    f.write('\n')
    print(len(Coords[0]))
    print(Atoms)
    # Now its when the party starts jeje
    
    for i in range(Atoms):
        #FOr each particle i have to print a special line
        id=i+1
        mol= int( i/Fil_L)+1
        if(i%Fil_L==0):
            atom_type=1
        elif(i%Fil_L==1):
            atom_type=2
        else:
            atom_type=3
        
        x_coord=Coords[0][i]
        y_coord=Coords[1][i]
        z_coord=0
        #Now i can write the line
        f.write('{} {} {} {} {} {}\n'.format(id,mol,atom_type,x_coord,y_coord,z_coord))
    #Nice you just wrote the first 
    
    
    
    f.write('\n')
    f.write('Bonds\n')
    f.write('\n')
    bonds_per_molecule=Fil_L-1
    for i in range(Bonds):
        id=i+1
        bond_type=1
        #First, how many bonds does a molecule have 
        #Bond number 5 is in the th molecule
        mol_num= int(i/bonds_per_molecule)
        
        part1= mol_num*Fil_L+ i%bonds_per_molecule +1
        part2=part1+1
        f.write('{} {} {} {}\n'.format(id,bond_type,part1,part2))



    f.write('\n')
    f.write('Angles\n')
    f.write('\n')

    #Now lets do the angles. for this one we will have a similar display i think
    #     
    angles_per_molecule=int(Fil_L-2)
    for i in range(Angles):
        id=i+1
        angle_type = 1
        mol_num = int(i/angles_per_molecule)
        part1= mol_num*Fil_L + i%angles_per_molecule +1
        part2= part1+1
        part3= part2 +1
        f.write('{} {} {} {} {}\n'.format(id,angle_type,part1,part2,part3))
        #We need to know the n
        
    f.close()




    return Config_dir


def main(Atoms, Fil_L,Box,f_swim,f_chem,f_chem2,sigma2,FField,Nsim):
    Config_dir='Configs/N_{}_L_{}_Lx_{}_Ly_{}_Nsim_{}_config_opt.txt'.format(Atoms,Fil_L,Box[1],Box[3],Nsim)
    g = open('./Configs/Logs/N_{}_L_{}_Lx_{}_Ly_{}_Nsim_{}_config_log_opt.txt'.format(Atoms,Fil_L,Box[1],Box[3],Nsim), 'w')
    start= time.time()
    #The input goes to the console and to the log
    String='Creating a simulation box with {} atoms in a box with limits x=({},{}), y=({},{})\n'.format(Atoms,Box[0],Box[1],Box[2],Box[3])
    print(String)
    g.write(String)

    n_density=Atoms/((Box[1]-Box[0])*(Box[3]-Box[2]))
    String='The number density is {}\n'.format(n_density )
    print(String )
    g.write(String)
    String='The area density is {}\n'.format(n_density*np.pi/4)
    print(String)
    g.write(String)
    # Atoms=int(300)
    # Fil_L=int(20)
    Bonds=int( (Fil_L-1)*Atoms/Fil_L)
    Angles=int( (Fil_L-2)*Atoms/Fil_L)
    Angle_types=1
    Atom_types=3
    Bond_types=1
    Config_file(Atoms,Fil_L,Box,Nsim,Bonds,Angles,Angle_types,Atom_types,Bond_types,Config_dir,g)
    Input_file(Atoms,Fil_L,Box,Config_dir,f_swim,f_chem,f_chem2,sigma2,FField,Nsim)
    end = time.time()
    elapsed=end-start
    String='This sample took {} seconds'.format(elapsed )
    print(String)
    g.write(String)    
    g.close()



# for vel in range(10):
#     main(Atoms,Fil_L,Box,10*vel+10,Nsim_ini)

f_chems1=[0.5]
f_chems2=[20,30]
FFields=[10]
# FFields=[0]
# NAtoms=[3000,6000,12000]
Nsimss=[20,21]
NAtoms=[15000]
sigma2s=[225]

for f_chem1 in f_chems1:
    for f_chem2 in f_chems2:
        for Atoms in NAtoms:
           for sigma2 in sigma2s:
                for Nsim in Nsimss:
                    main(Atoms,Fil_L,Box,f_swim,f_chem1,f_chem2,sigma2,FField,Nsim)

# for Nsim in range(Nsim_ini,Nsim_fin):
#     main(Atoms,Fil_L,Box,f_swim,f_chem,f_chem2,sigma2,FField,Nsim)


