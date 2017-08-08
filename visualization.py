import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.pyplot as plt
from numpy import linspace, meshgrid, sqrt, outer, sin, cos, ones, zeros, pi
from matplotlib.patches import Circle





def prepare_plot_space(x_lim,y_lim,z_lim):
    fig=plt.figure()
    ax = Axes3D(fig, azim=30, elev=30)
    ax.set_xlabel('x-axis')
    ax.set_xlim(x_lim)
    ax.set_ylabel('y-axis')
    ax.set_ylim(y_lim)
    ax.set_zlabel('z-axis')
    ax.set_zlim(z_lim)
    return(ax)


def plot_3D_cylinder(ax,radius, height, elevation=0, resolution=100, color='r', x_center = 0, y_center = 0):
    x =  linspace(x_center-radius, x_center+radius, resolution)
    z =  linspace(elevation, elevation+height, resolution)
    X, Z =  meshgrid(x, z)

    Y =  sqrt(radius**2 - (X - x_center)**2) + y_center # Pythagorean theorem

    ax.plot_surface(X, Y, Z, linewidth=0, color=color)
    ax.plot_surface(X, (2*y_center-Y), Z, linewidth=0, color=color)

    floor = Circle((x_center, y_center), radius, color=color)
    ax.add_patch(floor)
    art3d.pathpatch_2d_to_3d(floor, z=elevation, zdir="z")

    ceiling = Circle((x_center, y_center), radius, color=color)
    ax.add_patch(ceiling)
    art3d.pathpatch_2d_to_3d(ceiling, z=elevation+height, zdir="z")

def plot_sphere(ax,radius,elevation,color):

    u =  linspace(0, 2 *  pi, 100)
    v =  linspace(0,  pi, 100)
    x = radius *  outer( cos(u),  sin(v))
    y = radius *  outer( sin(u),  sin(v))+elevation
    z = radius *  outer( ones( size(u)),  cos(v))

    # Plot the surface
    ax.plot_surface(x, y, z, color=color)
    plt.show()
    #plt.close()




def plot_tree(ax,x,y,DBH,H,Crown_Length,Crown_Diameter):
    #plot stem
    plot_3D_cylinder(ax=ax,radius=DBH, height=H-Crown_Length, elevation=0, resolution=200, color="brown", x_center=x, y_center=y)
    #plot crown
    plot_3D_cylinder(ax=ax,radius=Crown_Diameter, height=Crown_Length, elevation=H-Crown_Length, resolution=100, color="green", x_center=x, y_center=y)



"""
plot_tree(ax=ax,x=30,y=20, DBH=1.5, H=15,Crown_Length=8,Crown_Diameter=5)
plot_tree(ax=ax,x=50,y=30, DBH=3, H=30,Crown_Length=12,Crown_Diameter=12)
plot_tree(ax=ax,x=80,y=80, DBH=1, H=10,Crown_Length=4,Crown_Diameter=5)
plt.show()
"""

def plot_NEE(s):
    #s=s*(-1)
    fig=plt.figure()
    ax1=plt.subplot2grid((1,1),(0,0))
    ax1.plot(range(len(s)),s,'-',label="NEE")
    ax1.plot([],[],linewidth=5,label="sink",color="b",alpha=0.5)
    ax1.plot([],[],linewidth=5,label="source",color="r",alpha=0.5)
    ax1.fill_between(range(len(s)),s,0,where=(s> zeros(len(s))),facecolor="b",alpha=0.5)
    ax1.fill_between(range(len(s)),s,0,where=(s< zeros(len(s))),facecolor="r",alpha=0.5)
    plt.title("Net Ecosystem Exchange")
    plt.xlabel("Time[yrs]")
    plt.ylabel("[to/ha/yr]")
    ax1.legend(loc='upper left', bbox_to_anchor=(1, 0.5))
    return fig
    #plt.savefig("NEE.png")
    #plt.close(fig)
    #plt.show()

def plot_emissions(dic):
    FastSoil=dic['soil_fast']
    SlowSoil=[f+s for f,s in zip(FastSoil,dic['soil_slow'])]
    DeadWood=[s+w for s,w in zip(SlowSoil,dic['dead_wood'])]
    LivingTrees=[d+l for d,l in zip(DeadWood,dic['living_trees'])]

    fig=plt.figure()
    ax1=plt.subplot2grid((1,1),(0,0))
    #ax1.plot(range(len(s)),s,'-',label="NEE")
    ax1.plot([],[],linewidth=5,label="Living Trees",color="green",alpha=0.5)
    ax1.plot([],[],linewidth=5,label="Dead Wood",color="orange",alpha=0.5)
    ax1.plot([],[],linewidth=5,label="Slow Soil",color="brown",alpha=0.5)
    ax1.plot([],[],linewidth=5,label="Fast Soil",color="yellow",alpha=0.5)

    length=len(dic['living_trees'])
    ax1.fill_between(range(length),LivingTrees,0,facecolor="green",alpha=0.8)
    ax1.fill_between(range(length),DeadWood,0,facecolor="orange",alpha=0.8)
    ax1.fill_between(range(length),SlowSoil,0,facecolor="brown",alpha=0.8)
    ax1.fill_between(range(length),FastSoil,0,facecolor="yellow",alpha=0.8)
    plt.title("Carbon Emmissions")
    plt.xlabel("Time[yrs]")
    plt.ylabel("Carbon [t]")
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=4)
    return fig

def plot_respiration(dic):
    Ft_col={"Cgpp":"forestgreen","Cr":"orange"}

    fig=plt.figure()
    ax1=plt.subplot2grid((1,1),(0,0))
    ax1.plot(dic["Cgpp"],linewidth=3,label="Cgpp",color=Ft_col["Cgpp"])
    ax1.plot(dic["Cr"],linewidth=3,label="Cr",color=Ft_col["Cr"])
    #plt.ylim(ymin=-10)
    plt.title("Respiration and Primary Production")
    plt.xlabel("Time[yrs]")
    plt.ylabel("Carbon")
    #plt.axhline(y=0)
    plt.legend(loc="upper right")
    return fig

def plot_populations(Pop):
    Ft_col={"FT1":"forestgreen","FT2":"orange","FT3":"saddlebrown","FT4":"mediumorchid","FT5":"crimson","FT6":"steelblue" }

    fig=plt.figure()
    ax1=plt.subplot2grid((1,1),(0,0))
    ax1.plot(Pop["FT1"],linewidth=3,label="FT1",color=Ft_col["FT1"])
    ax1.plot(Pop["FT2"],linewidth=3,label="FT2",color=Ft_col["FT2"])
    ax1.plot(Pop["FT3"],linewidth=3,label="FT3",color=Ft_col["FT3"])
    ax1.plot(Pop["FT4"],linewidth=3,label="FT4",color=Ft_col["FT4"])
    ax1.plot(Pop["FT5"],linewidth=3,label="FT5",color=Ft_col["FT5"])
    ax1.plot(Pop["FT6"],linewidth=3,label="FT6",color=Ft_col["FT6"])
    #plt.ylim(ymin=-10)
    plt.title("Functional Type Population")
    plt.xlabel("Time[yrs]")
    plt.ylabel("Population Size")
    #plt.axhline(y=0)
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=6)
    return fig
    #plt.show()
    #plt.close(fig)

def plot_stocks(Stocks):
    Soil=[s+f for s,f in zip(Stocks["Sslow"],Stocks["Sfast"])]
    DeadWood=[s+w for s,w in zip(Soil,Stocks["Dwood"])]
    AGB=[d+a for d,a in zip(DeadWood,Stocks["AGB"])]
    fig=plt.figure()
    ax1=plt.subplot2grid((1,1),(0,0))
    #ax1.plot(range(len(s)),s,'-',label="NEE")
    ax1.plot([],[],linewidth=5,label="AGB",color="green",alpha=0.5)
    ax1.plot([],[],linewidth=5,label="Dead Wood",color="orange",alpha=0.5)
    ax1.plot([],[],linewidth=5,label="Soil",color="brown",alpha=0.5)
    ax1.fill_between(range(len(AGB)),AGB,0,facecolor="green",alpha=0.5)
    ax1.fill_between(range(len(DeadWood)),DeadWood,0,facecolor="orange",alpha=0.5)
    ax1.fill_between(range(len(Soil)),Soil,0,facecolor="brown",alpha=0.5)
    plt.title("Carbon Stocks")
    plt.xlabel("Time[yrs]")
    plt.ylabel("Carbon [t]")
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=3)
    return fig
    #plt.show()
    #plt.close(fig)






def plot_ground_light(topology):
    plt.imshow(topology.surface[topology.dim_names["GroundLight"]], cmap=plt.get_cmap("gray"),vmin=0,vmax=860)
    plt.colorbar()
    plt.show()

def plot_tree_pos(Tree,output_dir,length,xlim,ylim,step):
    x=[]
    y=[]
    s=[]
    tp=[]
    Ft_col={"FT1":"forestgreen","FT2":"orange","FT3":"saddlebrown","FT4":"mediumorchid","FT5":"crimson","FT6":"steelblue" }

    for t in Tree.Instances.values():
        x.append(t.position[0])
        y.append(t.position[1])
        s.append(t.CA)
        tp.append(Ft_col[t.Ftype])

    fig=plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    #ax.autoscale(False)
    ax.set_xlim(-1,xlim+1)
    ax.set_ylim(-1,ylim+1)
    width, height = fig.canvas.get_width_height()
    p=length/width

    areas=[a/p for a in s]


    ax.scatter(x,y,s=areas,c=tp,alpha=0.5)
    plt.title("Simulation Year:"+str(step))
    plt.savefig(output_dir+"step: "+str(step)+".png")
    plt.close(fig)

def add_time_guide(figname,length,plot_func,**kwargs):
    """
    Plot the result of [plot_func] [length] times, with a vetical line moving along the x (time) axis.

    plot_func() must return a pylab fig.
    """
    for i in range(length):
        fig=plot_func(**kwargs)
        plt.axvline(x=i)
        plt.savefig(figname+str(i)+".png")
        plt.close(fig)
