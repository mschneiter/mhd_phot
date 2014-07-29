from huachoUtils import *
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import LogNorm


mpl.rc('xtick',labelsize=8)
mpl.rc('ytick',labelsize=8)
mpl.rc('text',usetex=True)
#mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})

def add_inner_title(ax, title, loc, size=None, **kwargs):
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.patheffects import withStroke
    if size is None:
        size = dict(size=plt.rcParams['ytick.labelsize'])
    at = AnchoredText(title, loc=loc, prop=size,
                      pad=0., borderpad=0.5,
                      frameon=False, **kwargs)
    ax.add_artist(at)
    at.txt._text.set_path_effects([withStroke(foreground="w", linewidth=2)])
    return at


nxtot=400
nytot=400
nztot=4
mpiX=4
mpiY=4
mpiZ=1
neqs=8


cv=1.5
gam=(cv+1.)/cv
vsc=1.
psc=1.
Bsc= np.sqrt(4.*math.acos(-1.))
firstrun=True
plt.ion()

path='/Users/esquivel/Desktop/Storage-Diable/MHD-EXO/OT-HLLD/BIN/'

x_ax=np.linspace(0,1,nxtot)

for nout in range(10,11):
    #  this is the ionization rate
    if (firstrun):
        denT= coplot3d(3,nztot/2,0,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)
#        vx  = coplot3d(3,nztot/2,1,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)/denT
#        vy  = coplot3d(3,nztot/2,2,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)/denT
#        vz  = coplot3d(3,nztot/2,3,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)/denT
#        Etot= coplot3d(3,nztot/2,4,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)
#        bx=   coplot3d(3,nztot/2,5,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)/Bsc
#        by=   coplot3d(3,nztot/2,6,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)/Bsc
#        bz=   coplot3d(3,nztot/2,7,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)/Bsc
        
#        Ekin= 0.5*(denT*(vx*vx+vy*vy+vz*vz))
#        Emag= 0.5*(bx*bx+by*by+bz*bz)
#        Pgas=  (Etot-Ekin-Emag)/cv
#        Temp= (Pgas)/(cv*denT)


    plt.figure(1)
    plt.clf()
    plt.imshow(denT,cmap='jet',origin='lower',vmin=0.05,vmax=0.5)
    plt.colorbar()
    plt.show()
    plt.savefig('HLLD/dens-OT-'+str(nout).zfill(3)+'.png',dpi=100,transparent=True,bbox_inches='tight')
    

#    plt.figure(2)
#    plt.clf()
#    plt.imshow(Ekin,cmap='jet',origin='lower' )
#    plt.colorbar()
#    plt.show()  

#    plt.figure(3)
#    plt.clf()
#    plt.imshow(Temp,cmap='jet',origin='lower')#,norm=LogNorm() )
#    plt.colorbar()
#    plt.show()  

        
#plt.ioff()

