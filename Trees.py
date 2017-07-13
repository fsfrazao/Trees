from Py_IBM import *
from plot_trees import plot_tree_pos
import pdb
from tables import *
import json

class Tree_World(World):
    """ Provide observer methods and control the schedule and execution of simulations.


    Args:
        topology (Tree_Grid instance): grid object on which agents will live.
        AET (int): actual evapotranspiration in the previous year in mm.
        Assumed to be constant every year.

    Attributes:
        Smort (float): Carbon of all trees that died in one year [tC]
        Sdead (float): Dead wood carbon stock [tC]
        tsdead (float): Anual decomposition rate
        tsdead_A (float): Proportion of decomposed dead wood that is released to the atmosphere.
        Sslow (float): Slow decomposition portion of soil stock [tC]
        tsdead_Sslow (float): Proportion of decomposed dead wood that is tranferred to the slow decomposition Soil stock
        self.tsslow_A (float):Proportion of slow decomposition soil stock that is released to the atmosphere.
        Sfast (float):Fast decomposition portion of soil stock [tC]
        tsfast_A (float):Proportion of fast decomposition soil stock that is released to the atmosphere.
        tsdead_Sfast (float): Proportion of decomposed dead wood that is tranferred to the fast decomposition Soil stock.
        Cgpp (float): Carbon captured in the gross primary productivity of the living forest
        Cr (float):Carbon released by the total respiration of the living forest (i.e. for maintenance and growth)

    Usage:
        >>>grid01=Tree_Grid(x_max=5,y_max=5,delta_h=0.5,patch_area=400,I0=860,lday=12,phi_act=360,k=0.7)
        >>>world01=Tree_World(topology=grid01)
        >>> print(world01)
        World object
    """

    def __init__(self, topology, AET=1300):
        super().__init__(topology)
        self.mortality_factor=1
        self.AET=AET
        self.Smort=0
        self.Sdead={"current":0,"previous":0}
        self.tsdead=min(1.0,np.power(10,(-1.4553+0.0014175*self.AET))/12)
        self.recruit_factor=1.0
        self.dwood_emission_p=0.7
        self.tsdead_A=self.dwood_emission_p*self.tsdead
        self.Sslow={"current":0,"previous":0}
        self.tsdead_Sslow=0.015*0.3*self.tsdead
        self.tsslow_A=1/750
        self.Sfast={"current":0,"previous":0}
        self.tsfast_A=1/15
        self.tsdead_Sfast=0.985*0.3*self.tsdead
        self.smort_log=[]
        self.Cgpp=0
        self.Cr=0
        self.NEE=[]
        self.Stocks={"Sslow":[0,],
                    "Sfast":[0,],
                    "Dwood":[0,],
                    "AGB":[0,],}
        self.tree_subclasses={"FT1":Tree_FT1,
                        "FT2":Tree_FT2,
                        "FT3":Tree_FT3,
                        "FT4":Tree_FT4,
                        "FT5":Tree_FT5,
                        "FT6":Tree_FT6}



    def run_simulation(self,n,logging_years,min_dbh,max_dbh,vol,h5file=None,plot_func=None,**kwargs):
        """Run the code in run_schedule() 'n' times.

        Args:
            n (int): Number of simulations to be run.
            logging_years (range): The years in which selective logging occurs.
            min_dbh (float): minimum dbh for a tree to be considered suitable for logging.
            max_dbh (float): maximum dbh for a tree to be considered suitable for logging.
            vol (float): total volume of timber to be extracted by logging event.
            h5file (str): path to the HDF5 file in which model outputs will be saved
            plot_func (function): a function that will be run every step. Usually a plotting function that saves a image file to the disk. Must accept a 'step' argument,which will receive the step number and can be used to generate the file name.
            **kwargs: List of keyword arguments to be passed to plot_func in addition to 'step'.

        Returns:
            None
        """
        if h5file:

            database=open_file(h5file,"r+")

            Pop_table=database.root.sim_1.sys_lvl.Pop
            Stocks_table=database.root.sim_1.sys_lvl.Stocks
            NEE_table=database.root.sim_1.sys_lvl.NEE
            #Emissions_table=database.root.sim_1.sys_lvl.Emissions


            pop_r=Pop_table.row
            stocks_r=Stocks_table.row
            nee_r=NEE_table.row
            #emissions_r=Emissions_table.row

        self.topology.update()
        for t in Tree.Instances.values():
            t.update()
        self.stablish_new_trees()

        for i in range(n):

            print("step=",i,"N=",len(Tree.Instances),"DBH>10",Tree.TreesAboveDBH(10.0), "DTrees", sum(Tree.DeadTrees.values()))
            if h5file:
               pop_r['FT1']=len(Tree_FT1.Indices)
               pop_r['FT2']=len(Tree_FT2.Indices)
               pop_r['FT3']=len(Tree_FT3.Indices)
               pop_r['FT4']=len(Tree_FT4.Indices)
               pop_r['FT5']=len(Tree_FT5.Indices)
               pop_r['FT6']=len(Tree_FT6.Indices)
               pop_r.append()

               stocks_r['agb']=Tree.Total_AGB()*.44
               stocks_r['dead_wood']=self.Sdead["current"]
               stocks_r['soil_slow']=self.Sslow["current"]
               stocks_r['soil_fast']=self.Sfast["current"]
               stocks_r.append()



            self.topology.update()
            for t in Tree.Instances.values():
                t.update()
            self.topology.update()


            if h5file:

               nee_r['dead_wood']=self.tsdead_A*self.Sdead["current"]
               nee_r['soil_slow']=self.tsslow_A*self.Sslow["current"]
               nee_r['soil_fast']=self.tsfast_A*self.Sfast["current"]
               nee_r['gpp']=self.Cgpp
               nee_r['living_trees']=self.Cr
               nee_r['smort']=self.Smort
               nee_r['t_deadwood_sslow']=self.tsdead_Sslow*self.Sdead["previous"]
               nee_r['t_sslow_sfast']=self.tsdead_Sfast*self.Sdead["previous"]



               nee=self.run_schedule(step=i,n=n,h5_database=database)

               nee_r['total_emissions']=-nee+self.Cgpp

               nee_r['nee']=nee
               nee_r.append()
               NEE_table.flush()
               Pop_table.flush()
               Stocks_table.flush()


            else:
                self.NEE.append(self.run_schedule(step=i,n=n))#,h5_database=database)
                self.Stocks["Sslow"].append(self.Sslow["current"])
                self.Stocks["Sfast"].append(self.Sfast["current"])
                self.Stocks["Dwood"].append(self.Sdead["current"])
                self.Stocks["AGB"].append(Tree.Total_AGB()*.44)


            if not plot_func==None: plot_func(step=i,**kwargs)

            if i in logging_years:
                ps=list(self.topology.trees_per_patch.keys())
                st=self.suitable_trees(ps,min_dbh,max_dbh)
                self.log_trees(st,vol)


            self.stablish_new_trees()
        if h5file: database.close()

    def run_schedule(self,step,n, h5_database=None):
        """Schedule of actions to be run in every time step

        Args:
            step (int): The current step. Received from run_simulation().
            h5_database (str): path to the HDF5 file in which model outputs will be saved

        Returns:
            (float) the Net Ecosystem Exchange (as calculater by 'calculate_NEE') for thar step.
        """
        #self.topology.update()

        list_of_patches=self.topology.trees_per_patch.keys()
        cohorts=self.topology.cohorts(list_of_patches)
        FTs=list(cohorts.keys())


        for ft in FTs:
            ages=list(cohorts[ft].keys())
            for age in ages:
                base_tree_index=cohorts[ft][age][0]
                base_tree=Tree.Instances.get(base_tree_index)
                base_tree.increase_DBH()
                base_tree.age+=1
                attributes=base_tree.update()

                for tree_index in cohorts[ft][age][1:]:
                    tree=Tree.Instances.get(tree_index)
                    tree.import_attributes(attributes)



        #print("DBH:",Tree.Instances[125].DBH,"AGB:",Tree.Instances[125].AGB )

        #trees=[k for k in Tree.Instances.keys()]
        #for i in trees:
        #    t=Tree.Instances.get(i)
            #print("dbh:",t.id)
        #    t.increase_DBH()

        #

        if h5_database:
            Ind_table=h5_database.root.sim_1.ind_lvl.Ind
            ind_r=Ind_table.row
            trees=[k for k in Tree.Instances.keys()]
            for i in trees:
               t=Tree.Instances.get(i)
            #    #print("up:",t.id)
            #    t.update()
            #    t.age+=1


               ind_r["step"]=step
               ind_r["ind_id"]=t.id
               ind_r["pos_x"]=round(t.position[0],2)
               ind_r["pos_y"]=round(t.position[1],2)
               ind_r["age"]=t.age
               ind_r["dbh"]=t.DBH
               ind_r["agb"]=t.AGB
               ind_r["height"]=t.H
               ind_r["crown_area"]=t.CA
               ind_r["FT"]=t.Ftype
               ind_r.append()

            Ind_table.flush()

        self.Cgpp=Tree.Cgpp()
        self.Cr=Tree.Cr()
        Tree.Mortality()
        #Tree.BackgroundMortality()
        # Tree.CrowdingMortality()
        # Tree.DamageMortality()
        return self.calculate_NEE()


    def create_agents(self,FT,n,pos=None, ids=None, DBHs=None,ages=None, **kwargs):
        """ Creates 'n' trees of type 'FT'. A list of positions, ids and DBHs can be provided (Useful to recreate trees from a previous model run, for example).

        Note: overwrites methods create_agents from class Agent.

        Args:
            agent_type (class): the type (class) of agents to be created
            n (int): number of agents to be created
            pos (list,None): a list of 'n' tuples '(x,y)' with the coordinates in which agents will be created. If set to None, agents are created in random positions.
            ids (list): a list of 'n' ints to be used as tree ids.
            DBHs: a list of 'n' floats to be used as tree DBHs.
            kwargs: the arguments to be passed to the 'agent_type' class. All agents will be created with this same set of arguments.

        Returns:
            None.
        """
        if pos==None:
            pos=[(np.random.rand()*self.topology.x_max,
            np.random.rand()*self.topology.y_max)
            for i in range(n)]

        if ids==None:
            for i in range(n):
                FT(position=pos[i],**kwargs)

        else:
            for i in range(n):
                FT(position=pos[i],dbh=DBHs[i], id=ids[i],age=ages[i],**kwargs)


    def seedbank_from_file(self, input_file):
        with open(input_file,'r') as json_file:
            self.topology.seedbank=json.load(json_file)



    def seedbank_to_file(self,output_file):
        with open(output_file,'w') as json_file:
            json.dump(self.topology.seedbank,json_file)


    def model_status_from_file(self, input_file):
        with open(input_file,'r') as json_file:
            trees=json.load(json_file)

        for ft in self.tree_subclasses.keys():
            ages=[]
            positions=[]
            ids=[]
            DBHs=[]
            for tree in trees[ft]:
                positions.append(tuple(tree['position']))
                ids.append(tree['id'])
                DBHs.append(tree['DBH'])
                ages.append(tree['age'])

            self.create_agents(self.tree_subclasses[ft],len(trees[ft]),pos=positions,DBHs=DBHs,ids=ids, ages=ages, world=self)

        self.topology.update()
        for t in Tree.Instances.values():
            t.update()
            t.AGB=t.AGB_from_DBH(t.DBH/100)
        Tree.ID=max(Tree.Instances.keys())+1


    def model_status_to_file(self,output_file):
        trees_per_type={k:[] for k in self.tree_subclasses.keys()}
        for t_id,t in Tree.Instances.items():
            trees_per_type[t.Ftype].append({"id":t_id, "position":t.position,"DBH":t.DBH,"age":t.age})

        with open(output_file,'w') as json_file:
            json.dump(trees_per_type,json_file)


    def suitable_trees(self, patches, min_dbh, max_dbh):
        """Search for suitable trees for logging in the given patches.

        Args:
            patches (list): A list of patches (tuples in the format (x,y))
            min_dbh (float): minimum dbh for a tree to be considered suitable for logging.
            max_dbh (float): max dbh for a tree to be considered suitable for logging.

        Returns:
            A dictionary with patches as keys ((x,y)) and lists of suitable tree ids as values
        """

        suitables_per_patch={}
        for p in patches:
            suitables_per_patch[p]=[t for t  in self.topology.trees_per_patch[p] if Tree.Instances[t].DBH>=min_dbh and Tree.Instances[t].DBH<=max_dbh]

        suitables_per_patch={k:v for k,v in suitables_per_patch.items() if v}
        return suitables_per_patch

    def log_trees(self, suitable_trees, total_vol):
        """Randomly log 'suitable_trees' until the 'total_vol' is reached.

        Args:
            suitable_trees (dict): A dictionary with patches as keys ((x,y)) and lists of suitable tree ids as values (resulting from suitable_trees())
            total_vol (float): the total vol of timber to be logged

        Returns:
            None
        """

        total_logged=0
        while total_logged<total_vol and len(suitable_trees):
            ps=list(suitable_trees.keys())
            p=ps[np.random.randint(len(ps))]
            vols={t:Tree.Instances[t].calculate_volume() for t in suitable_trees[p]}
            #pdb.set_trace()
            biggest_tree=max(vols, key=lambda x: vols[x])
            total_logged+=vols[biggest_tree]
            Tree.Instances[biggest_tree].log_me()
            suitable_trees[p].remove(biggest_tree)
            if not suitable_trees[p]:
                del suitable_trees[p]



    def define_seedbank(self,Nseed):
        """Randomly determines how many seed each patch will receive.

        Args:
            Nseed (int):The total number of seed to be distributed among all patches.

        Returns:
            A dictionary in which keys are patches ('(x,y)') and values are the number of seeds
        """

        seedbank={}
        for s in range(Nseed):
            p=(np.random.randint(self.topology.x_max),np.random.randint(self.topology.y_max))
            if seedbank.get(p)==None:
                seedbank[p]=1
            else:
                seedbank[p]+=1

        return seedbank

    def seeds_pos(self,est_seeds):
        """Set a random position within the patch for each seed in est_seeds.

        Args:
            est_seed (dict):A dictionary with patches as keys ('(x,y)') and the number of established seeds as values.

        Returns:
            A list of seed positions within patches.
        """

        pos=[]
        for i in est_seeds.items():
            for j in range(i[1]):
                seed_pos=i[0]+np.random.rand(2)
                pos.append(self.topology.in_bounds(seed_pos))
        return pos



    def germinate_suitable_seeds(self,FT,pos):
        self.create_agents(self.tree_subclasses[FT],len(pos),pos=pos, world=self, dbh=2)


    def stablish_new_trees(self):
        """Create new tree individuals. Use 'define_seedbank()', 'topology.seed_establishment' and 'seeds_pos'.

            Note: This method can be substitued by an equivalent in order to incorporate more realistic representations of seed dispersal.

        Returns:
            None.
        """


        seedbank_1=self.define_seedbank(Nseed=int(Tree_FT1.Nseed*self.topology.total_area))
        est_seeds_1=self.topology.seed_establishment(seedbank_1,lmin=Tree_FT1.Iseed,h0=Tree_FT1.h0,h1=Tree_FT1.h1)
        pos=self.seeds_pos(est_seeds_1)
        self.topology.seedbank["FT1"]=pos
        self.germinate_suitable_seeds(FT="FT1",pos=pos)
        #self.create_agents(Tree_FT1,n,pos=pos, world=self,dbh=2)

        seedbank_1=self.define_seedbank(Nseed=int(Tree_FT2.Nseed*self.topology.total_area))
        est_seeds_1=self.topology.seed_establishment(seedbank_1,lmin=Tree_FT2.Iseed,h0=Tree_FT2.h0,h1=Tree_FT2.h1)
        pos=self.seeds_pos(est_seeds_1)
        self.topology.seedbank["FT2"]=pos
        self.germinate_suitable_seeds(FT="FT2",pos=pos)

        seedbank_1=self.define_seedbank(Nseed=int(Tree_FT3.Nseed*self.topology.total_area))
        est_seeds_1=self.topology.seed_establishment(seedbank_1,lmin=Tree_FT3.Iseed,h0=Tree_FT3.h0,h1=Tree_FT3.h1)
        pos=self.seeds_pos(est_seeds_1)
        self.topology.seedbank["FT3"]=pos
        self.germinate_suitable_seeds(FT="FT3",pos=pos)

        seedbank_1=self.define_seedbank(Nseed=int(Tree_FT4.Nseed*self.topology.total_area))
        est_seeds_1=self.topology.seed_establishment(seedbank_1,lmin=Tree_FT4.Iseed,h0=Tree_FT4.h0,h1=Tree_FT4.h1)
        pos=self.seeds_pos(est_seeds_1)
        self.topology.seedbank["FT4"]=pos
        self.germinate_suitable_seeds(FT="FT4",pos=pos)

        seedbank_1=self.define_seedbank(Nseed=int(Tree_FT5.Nseed*self.topology.total_area))
        est_seeds_1=self.topology.seed_establishment(seedbank_1,lmin=Tree_FT5.Iseed,h0=Tree_FT5.h0,h1=Tree_FT5.h1)
        pos=self.seeds_pos(est_seeds_1)
        self.topology.seedbank["FT5"]=pos
        self.germinate_suitable_seeds(FT="FT5",pos=pos)

        seedbank_1=self.define_seedbank(Nseed=int(Tree_FT6.Nseed*self.topology.total_area))
        est_seeds_1=self.topology.seed_establishment(seedbank_1,lmin=Tree_FT6.Iseed,h0=Tree_FT6.h0,h1=Tree_FT6.h1)
        pos=self.seeds_pos(est_seeds_1)
        self.topology.seedbank["FT6"]=pos
        self.germinate_suitable_seeds(FT="FT6",pos=pos)





    def calculate_Sdead(self):
        """ Calculate the amount of carbon (tC) in the dead wood stock based
            on the trees that died in the previous year (Smort).

            Returns:
                None.

        """
        #self.Smort/=2
        self.Sdead["current"]=max(0,self.Sdead["previous"]+(self.Smort)-self.tsdead*self.Sdead["previous"])
        self.Smort=0

    def calculate_Sslow(self):
        """ Calculate the amount of carbon (tC) in the slow decomposition soil
            stock.

            Calculation is based on decomposition rates for the dead wood (tsdead_Sslow) and the rate with which carbon is transferred from this stock to the atmosphere (tsslow_A).

            Returns:
                None.
        """

        self.Sslow["current"]=self.Sslow["previous"]+self.tsdead_Sslow*self.Sdead["previous"]-self.tsslow_A*self.Sslow["previous"]

    def calculate_Sfast(self):
        """ Calculate the amount of carbon (tC) in the fast decomposition soil
            stock.

            Calculation is based on decomposition rates for the dead wood (tsdead_Sfast) and the rate with which carbon is transferred from this stock to the atmosphere (tsfast_A).

            Returns:
                None.
        """

        self.Sfast["current"]=self.Sfast["previous"]+self.tsdead_Sfast*self.Sdead["previous"]-self.tsfast_A*self.Sfast["previous"]

    def calculate_NEE(self):
        """
        Return the Net Ecosystem Exchange (NEE) for the previous timestep.

        This is all the carbon that was absorbed by living trees (Cgpp) minus
        what was emitted to the atmosphere by the living trees (Cr) and
        the decomposition of the dead wood and soil stocks (tsdead_A+tsslow_A+tsfast_A).

        Returns:
            None.

        """

        #self.smort_log.append(self.Smort)

        self.calculate_Sdead()
        self.calculate_Sslow()
        self.calculate_Sfast()


        NEE= self.Cgpp-self.Cr-(self.tsdead_A*self.Sdead["previous"])-(self.tsslow_A*self.Sslow["previous"])-(self.tsfast_A*self.Sfast["previous"])

        self.Sdead["previous"]=self.Sdead["current"]
        self.Sslow["previous"]=self.Sslow["current"]
        self.Sfast["previous"]=self.Sfast["current"]



        return NEE




class Tree_Grid(Rectangular_Grid):
    """ Represent space as a grid.


    Args:
        x_max (int): number of horizontal cells
        y_max (int): number of vertical cells
        patch_area (float): are of each patch (cells) in squared meters
        I0 (float): Incoming irradiance on top of canopy [umol photon / m^2s]
        k (float): Light extinction coefficient.
        delta_h (float): height of each vertical layer.
        phi_act (float):photosynthetic active period.
        lday (int):mean day length during the vegetation period phi_act.
        ndim (int): number of dimensions. One per environmental variable.
        dim_names (list): a list of 'ndim' strings representing the names of dimensions.

    Attributes:
        LAIs (dict): dictionary with patches as keys and a list of averaged Leaf Area Indexes (one for each vertical layer) as values.
        CCAs (dict):dictionary with patches as keys a a list of Cummulative Crown Areas (one for each vertical layer) as values.
        trees_per_patch (list): dictionary with patches ('(x,y)') as keys and a list of tree ids as values
        total_area(float): area of all simulated patches in ha


    """

    def __init__(self,x_max,y_max,patch_area,I0,k,delta_h,phi_act,lday,ndim=5,dim_names=["GroundLight","LianaCoverage","Temperature","WindSpeed","Humidity"]):
        super().__init__(x_max,y_max,ndim=ndim,dim_names=dim_names)
        self.patch_area=patch_area
        self.delta_h=delta_h
        self.I0=I0
        self.k=k
        self.lday=lday
        self.phi_act=phi_act
        self.LAIs={}
        self.CCAs={}
        self.total_area=self.x_max * self.y_max *self.patch_area/10000
        #self.seedbanks={}
        self.trees_per_patch={(x,y):[] for x in range(0,self.x_max) for y in range(0,self.y_max)}
        self.FTs=["FT1","FT2","FT3","FT4","FT5","FT6"]
        self.seedbank={k:[] for k in self.FTs}





    def exeeding_CCA(self):
        patches=[patch for patch,layers in self.CCAs.items() if any(layer>self.patch_area for layer in layers )]

        return patches


    def update(self):
        """ Execute the methods that calculate values for each dimension of the grid.

        Returns:
            None.

        """

        self.LAIs=self.patch_based_LAI(Tree.TreesPerPatch())
        self.CCAs=self.patch_based_CCA(Tree.TreesPerPatch())
        self.ground_level_light()
        self.surface[self.dim_names["LianaCoverage"]]=self.calculate_liana_coverage()

    def seed_establishment(self,seedbank,lmin,h0,h1):
        """Check if the seeds in the seedbank have the light and space conditions to germinate.

            Note: All seeds in the given seedbank are assumed to have the same conditions. Multiple calls to this method passing different seedbanks allow the representation of different requirements.

            All established individuals are initiated with a predefined DBH (the same for all individuals), but initial height is calculated based on functional relationships and only if there is enough space in that layer a new individual will be stablished.

        Args:

            seedbank (dict): dictionary with patches ('(x,y)') and number of seed as values.
            lmin (float): minimum light intensity at ground level for seed to germinate.
            h0 (float): allometric paramenter used to calculate initial height
            h1 (float): allometric paramenter used to calculate initial height

        Returns:
            est_seeds: a dictionary with patches as keys ('(x,y)') and number of established seeds as values.
        """
        est_seeds={}
        for patch in seedbank.keys():
            light=self.check_light(lmin,patch)
            space=self.check_CCA(h0,h1,patch)
            if light and space:
                est_seeds[patch]=seedbank[patch]
            else:
                est_seeds[patch]=0
        return est_seeds



    def check_light(self,lmin, patch):
        """ Check if the light intensity at ground level meets the needs for seed germination.

            Note: Ligth intensity is homogeneized at the patch level.

        Args:
            lmin (float): minimum light intensity at ground for seeds to germinate.
            patch (tuple): (x,y), patch  that will have the light conditions assessed.

        Returns:
            Bool. True if light level is above minimum.
        """
        return self.surface[self.dim_names["GroundLight"]][patch]>lmin

    def check_CCA(self,h0,h1,patch):
        """ Check if the there is enough space forthe new plants to be established

            Note: Cummulative Crown Area is homogeneized at the patch level

        Args:
            h0 (float): allometric paramenter used to calculate initial height
            h1 (float): allometric paramenter used to calculate initial height
            patch (tuple): (x,y), patch  that will have the light conditions assessed.

        Returns:
            Bool. True if there is enough space.
        """
        Hmin=h0*2**h1
        l=int(np.floor(Hmin/self.delta_h))
        if self.CCAs.get(patch)==None:
            return True
        else:
            return self.CCAs[patch][l]<1.0


    def ground_level_light(self):
        """Calculate the amount of light that reaches the ground.

            The value per patch is stored in the "GroundLight" dimension.
            Values are homogeneized at the patch level.

            Returns:
                None.
        """
        self.surface[self.dim_names["GroundLight"]]=self.I0
        for patch in self.LAIs.keys():
            self.surface[self.dim_names["GroundLight"]][patch]=self.I0*np.exp(-self.k*sum(self.LAIs[patch]))



    def patch_based_LAI(self, trees_per_patch):
        """
        Calculate the Leaf Area Index for each patch in the grid.


        Args:
            trees_per_patch (dict): dictionary in which keys are patch coordinates (x,y) and values are  a list of the tree ids in that patch.

        Returns:
            L_per_layer(dict): dictionary in which the keys are patch coordinates('(x,y)')-the same as the ones in trees_per_patch- and values are a list of the LAI per vertical layer
        """

        L_per_layer={}
        nlayers=int(np.ceil(np.ceil(Tree.MaxHeight())/self.delta_h))
        for k in trees_per_patch.keys():
            L_per_layer[k]=[]
            lower_limit=int(min([Tree.Instances.get(t).lmin for t in trees_per_patch[k]]))
            upper_limit=int(max([Tree.Instances.get(t).lmax for t in trees_per_patch[k]]))


            for i in range(nlayers+1):

                L_per_layer[k].append((1/self.patch_area)*sum([Tree.Instances.get(tree_id).L_mean for tree_id in trees_per_patch[k] if lower_limit<=i and upper_limit>=i]))

                #upper and lower limits
                #Tree.Instances.get(tree_id).lmin<=i*self.delta_h and Tree.Instances.get(tree_id).lmax>=i*self.delta_h
        return L_per_layer

    def patch_based_CCA(self, trees_per_patch):
        """Calculate the Leaf Area Index for each patch in the grid.

        Args:
            trees_per_patch (dict): dictionary in which keys are patches ('(x,y)')' and values are lists of the trees ids in that patch.

        Returns:
            L_per_layer(dict): dictionary in which the keys are the same as the ones in trees_per_patch and values are a list of the CCA per vertical layer

            """
        L_per_layer={}
        nlayers=int(np.ceil(np.ceil(Tree.MaxHeight())/self.delta_h))
        for k in trees_per_patch.keys():
            L_per_layer[k]=[]
            lower_limit=int(min([Tree.Instances.get(t).lmin for t in trees_per_patch[k]]))
            upper_limit=int(max([Tree.Instances.get(t).lmax for t in trees_per_patch[k]]))

            for i in range(nlayers+1):

                L_per_layer[k].append((1/self.patch_area)*sum([Tree.Instances.get(tree_id).CA for tree_id in trees_per_patch[k] if lower_limit<=i and upper_limit>=i]))
        return L_per_layer

    def trees_in_patches(self,patches):
        """ Make a list of tree ids and their respective positions.

        Args:
            patch (list): a list of tuples representing the patch coordinates '(x,y)'.

        Returns:
            id_pos (list):a list tuples with tree ids and positions '(id,(x_pos,y_pos))'

        """


        #one list of ids for each patch
        ids=[self.trees_per_patch[p] for p in patches]
        #one single list of tuples with (id, pos)
        id_pos=[(tree,Tree.Instances.get(tree).position) for patch in ids for tree in patch]

        return id_pos


    def trees_per_age_groups(self, list_of_trees):
        """
        Return a dictionary in which keys are a age group and values are lists containing the id of the respective trees.

        list_of_trees is a list with the id of trees.
        """

        age_groups={}
        for t in list_of_trees:
            if age_groups.get(Tree.Instances[t].age)==None:
                age_groups[Tree.Instances[t].age]=[Tree.Instances[t].id]
            else:age_groups[Tree.Instances[t].age].append(Tree.Instances[t].id)

        return age_groups



    def trees_per_type(self, list_of_trees):
        """Organize trees by functional type.

        Args:
            list_of_trees (list): a list of tree ids.

        Returns:
            type_groups(dict): a dictionary in which keys are a Functional
            Types and values are lists containing the id of the respective trees.
        """

        type_groups={}
        for t in list_of_trees:
            if type_groups.get(Tree.Instances[t].Ftype)==None:
                type_groups[Tree.Instances[t].Ftype]=[Tree.Instances[t].id]
            else:type_groups[Tree.Instances[t].Ftype].append(Tree.Instances[t].id)

        return type_groups


    def cohorts(self, patches):
        """ Assemble cohorts based on tree type and age.

            Uses the trees_per_type() and trees_by_age_groups() methods to
            organize a dictionary of dictionaries containing the indices of
            trees in each age group within each type



        Args:
            Patches (list): a list of tuples representing patch coordinatas
            [(x1,y1),(x2,y2)...]

        Returns:
            cohorsts(dict): a dictionary with trees organized by type and age:
             {FunctionalType1:{AgeGroup1:[tree1,tree2,tree3...]}}
        """

        trees_in_patch=self.trees_in_patches(patches)
        list_of_trees=[i[0] for i in trees_in_patch]
        FT=self.trees_per_type(list_of_trees)
        cohorts={}
        for group in FT.keys():
          cohorts[group]=self.trees_per_age_groups(FT[group])

        return cohorts




    def calculate_liana_coverage(self):
        """ Calculate the coverage(%) os lianas in on patch based on the previous coverage, temperature, humidity and wind speed

        Returns:
                None
        """
        return np.zeros((self.x_max,self.y_max))

    def calculate_liana_biomass(self):
        """
        Calculate the biomass of lianas in one patch based on the coverage.
        """
        pass


class Tree(Agent):
    """ Represent individual trees.

    Args:
        position (tuple): (x,y) coordinates of where tree will be stablished.
        world (world object): world where tree will live.
        dbh (float): Diameter at Breast Height (in cm).
        h0 (float): Height-stem diameter relationship parameter.
        h1 (float): Height-stem diameter relationship parameter.
        cl0 (float): Crown Length-Height relationship parameter.
        cd0 (float): Crown Diameter-Stem Diameter relationship parameter.
        cd1 (float): Crown Diameter-Stem Diameter relationship parameter.
        cd2 (float): Crown Diameter-Stem Diameter relationship parameter.
        rho (float): wood density in t of Organic Dry Matter/cubic meter.
        sigma (float): ratio of total aboveground biomass to stem biomass.
        f0 (float): form factor-stem diameter relationship parameter.
        f1 (float): form factor-stem diameter relationship parameter.
        l0 (float): Leaf Area Index (LAI)-Stem relationship parameter.
        l1 (float): Leaf Area Index (LAI)-Stem relationship parameter.
        m (float): Light trasnmission coefficient.
        alpha (float): quantum eciency; initial slope of the type specic light response curve.
        pmax (float): maximum leaf gross photosynthetic rate.
        rg (float): fraction of gross primary production available for growth that is attributed to growth respiration.
        deltaDmax (float): maximum diameter increment.
        DdeltaDmax (float): % of Dmax that reaches deltaDmax.
        age (int): age (years).
        Ftype (str): functional type.
        Mb (float): backgroung mortality probability.
        Dfall (float): mimimum dbh for a tree to fall.
        pfall (float): probabbility that a tree with dbh>Dfall will fall.
        m_max (float): maximum size-dependent mortality for small trees.
        Dmort (float): DBH up to which mortality is increased (for small trees).

    Attributes:

        ID (int):
        Instances (dict):
        DeadTrees (dict):

        patch (tuple): (x,y) patch in which tree is located.
        f (float):
        k (float):
        lday (float):
        phi_act (float):

        H (float): Height (m).
        CL (float): Crown length (m).
        CD (float): Crown diameter (m)
        CA (float): Crown area (sq. m)
        AGB (float): Above ground biomass (tODM)
        LAI (float): Leaf area index.
        lmax (float):
        lmin (float):
        L_mean (float):
        Li (float):
        Mc (float):
        basic_m (float):
        Mb (float):
        Md (float):
        GPP (float): Gross primary production in this time step.
        Rm (float): Respiration maintenance in this timestep.
        position (tuple): (x,y) coordinates of where tree is located.
        world (world object): world where tree lives.
        dbh (float): Diameter at Breast Height (in cm).
        h0 (float): Height-stem diameter relationship parameter.
        h1 (float): Height-stem diameter relationship parameter.
        cl0 (float): Crown Length-Height relationship parameter.
        cd0 (float): Crown Diameter-Stem Diameter relationship parameter.
        cd1 (float): Crown Diameter-Stem Diameter relationship parameter.
        cd2 (float): Crown Diameter-Stem Diameter relationship parameter.
        rho (float): wood density in t of Organic Dry Matter/cubic meter.
        sigma (float): ratio of total aboveground biomass to stem biomass.
        f0 (float): form factor-stem diameter relationship parameter.
        f1 (float): form factor-stem diameter relationship parameter.
        l0 (float): Leaf Area Index (LAI)-Stem relationship parameter.
        l1 (float): Leaf Area Index (LAI)-Stem relationship parameter.
        m (float): Light trasnmission coefficient.
        alpha (float): quantum eciency; initial slope of the type-specic light response curve.
        pmax (float): maximum leaf gross photosynthetic rate.
        rg (float): fraction of gross primary production available for growth that is attributed to growth respiration.
        deltaDmax (float): maximum diameter increment.
        DdeltaDmax (float): % of Dmax that reaches deltaDmax.
        age (int): age (years).
        Ftype (str): functional type.
        Mb (float): backgroung mortality probability.
        max_height (float): maximum height this tree can reach.
        Dfall (float): mimimum dbh for a tree to fall.
        pfall (float): probabbility that a tree with dbh>Dfall will fall.
        m_max (float): maximum size-dependent mortality for small trees.
        Dmort (float): DBH up to which mortality is increased (for small trees).

    """


    ID=0
    Instances={}
    DeadTrees={"FT1":0,"FT2":0,"FT3":0,"FT4":0,"FT5":0,"FT6":0}

    @classmethod
    def TreesAboveDBH(self,dbh):
        """Calculate the total number of trees with DBH above the speciefied threshold.

        Returns:
            (int): Number of trees.
        """

        trees=[k for k in self.Instances.keys() if self.Instances.get(k).DBH>=dbh]
        return len(trees)

    @classmethod
    def Total_AGB(self):
        """ Calculate the total Above Ground Biomass for all living trees.

        Returns:
            (float): total AGB in t of Organic Dry matter (tODM).
        """

        Tagb=0
        trees=[k for k in self.Instances.keys()]
        for i in trees:
            t=Tree.Instances.get(i)
            Tagb+=t.AGB
        return Tagb

    @classmethod
    def Cgpp(self):
        """ Calculate the total Carbon absorbed as Gross Primary Production.

        Returns:
            (float): Total Carbon absorbed.
        """


        Cgpp=0
        trees=[k for k in self.Instances.keys()]
        for i in trees:
            t=Tree.Instances.get(i)
            Cgpp+=t.GPP
        return 0.44*Cgpp

    @classmethod
    def Cr(self):
        """ Calculate the total Carbon released through respiration.

        Returns:
            (float): Total Carbon released by all living trees.
        """
        Cr=0
        trees=[k for k in self.Instances.keys()]
        for i in trees:
            t=Tree.Instances.get(i)
            Cr+=(t.Rm+t.rg*(t.GPP-t.Rm))
        return 0.44*Cr

    @classmethod
    def Mortality(self):
        """
        Kill trees by calling their stochastic_M() method.
        """
        trees=[k for k in Tree.Instances.keys()]
        for t in trees:
             if Tree.Instances.get(t).stochastic_M(): Tree.Instances.get(t).remove_me()


    @classmethod
    def BackgroundMortality(self):
        """
        Kill trees by calling their stochastic_Bm() method.
        """
        trees=[k for k in self.Instances.keys()]
        for t in trees:
             if Tree.Instances.get(t).stochastic_Bm(): Tree.Instances.get(t).remove_me()

    @classmethod
    def CrowdingMortality(self):
         """
         Kill trees by calling their stochastic_Cm() method.
         """
         trees=[k for k in self.Instances.keys()]
         for t in trees:
              if Tree.Instances.get(t).stochastic_Cm(): Tree.Instances.get(t).remove_me()

    @classmethod
    def DamageMortality(self):
          """
          Kill trees by calling their stochastic_Dm() method.
          """
          trees=[k for k in self.Instances.keys()]
          for t in trees:
               if Tree.Instances.get(t).stochastic_Dm(): Tree.Instances.get(t).remove_me()


    @classmethod
    def ActiveAgents(self):
        return [ag.id for ag in self.Instances.values() if ag.active == True]

    @classmethod
    def TreesPerPatch(self):
        """Identify which trees are in each patch.

        Returns:
            (dict): keys are patch coordinates '(x,y)' and values are lists of trees ids in that patch.

        """
        patches={}
        for t in self.Instances.values():
            if patches.get(t.patch)==None:
                patches[t.patch]=[t.id]
            else:
                patches[t.patch].append(t.id)
        return patches

    @classmethod
    def MaxHeight(self):
        """ Finds the height of the tallest tree in the current time step.

        Returns:
            (float): height of the tallest tree.

        """
        return max([t.H for t in self.Instances.values()],default=0)


    def __repr__(self):
        return "Tree type: {0}".format(self.Ftype)

    def __init__(self,position,world,dbh,wood_density,h0,h1,cl0,cd0,cd1,cd2,rho,sigma,f0,f1,l0,l1,m,alpha,pmax,max_agb,rg,deltaDmax,DdeltaDmax,age,Ftype,Mb,Dfall=45, pfall=0.3,m_max=0.12,Dmort=10,max_dbh=None,id=None):

        super().__init__(position,world,id=id)
        #self.AddInstance(self.id)
        #Tree.Instances[Ftype+' '+str(self.id)]=self

        self.patch=(int(np.floor(self.position[0])),int(np.floor(self.position[1])))
        self.DBH=dbh
        self.age=age
        self.Ftype=Ftype
        self.wood_density=wood_density
        #type specific parameter used to calculate tree height
        self.h0=h0
        self.h1=h1
        #type specific parameter used to calculate crown length
        self.cl0=cl0
        #type specific parameter used to calculate crown diameter
        self.cd0=cd0
        self.cd1=cd1
        self.cd2=cd2
        #type specific parameter used to calculate aboveground biomass
        self.rho=rho
        self.sigma=sigma
        self.f0=f0
        self.f1=f1
        self.f=self.calculate_f()
        #type specific parameter used to calculate the leaf area index
        self.l0=l0
        self.l1=l1
        self.k=self.world.topology.k
        self.lday=self.world.topology.lday
        self.phi_act=self.world.topology.phi_act
        self.m=m
        self.alpha=alpha
        self.pmax=pmax
        self.m_max=m_max
        self.Dmort=Dmort
        self.Dfall=Dfall
        self.pfall=pfall
        self.Dmax=max_dbh/100.
        self.deltaDmax=deltaDmax
        self.DdeltaDmax=DdeltaDmax
        self.alpha0=self.calculate_alpha0()
        self.alpha1=self.calculate_alpha1()
        self.rg=rg

        self.H=0
        self.CL=0
        self.CD=0
        self.CA=0
        self.AGB=0
        self.LAI=0
        self.lmax=0
        self.lmin=0
        self.L_mean=0
        self.Li=0
        self.Mc=0
        self.basic_m=Mb
        self.Mb=0
        self.Md=0
        self.GPP=0 #Gpp in this time step
        self.Rm=0 #Rm in this timestep
        #self.update()
        self.world.topology.trees_per_patch[self.patch].append(self.id)


    def update(self):
        """ Updates geometry related parameters and mortality probabilities.
        """
        self.H=self.H_from_DBH()
        self.CL=self.CL_from_H()
        self.CD=self.CD_from_DBH()
        self.CA=self.CA_from_CD()
        if self.age==0:
            self.AGB=self.AGB_from_DBH(D=self.DBH/100)
        self.LAI=self.LAI_from_DBH()

        self.lmax=self.H/self.world.topology.delta_h #eq.31
        self.lmin=(self.H-self.CL)/self.world.topology.delta_h #eq.32
        self.L_mean=(self.LAI*self.CA)/(self.lmax-self.lmin if self.lmax-self.lmin!=0 else 0.5 )#eq.34
        self.Li=self.Lind()
        self.Mc=self.calculate_Rc()
        self.Mb=max(0,self.basic_m+self.calculate_Ms())
        self.Md=0


        attributes={"DBH":self.DBH,
                "H":self.H,
                "CL":self.CL,
                "CD":self.CD,
                "CA":self.CA,
                "AGB":self.AGB,
                "LAI":self.LAI,
                "lmax":self.lmax,
                "lmin":self.lmin,
                "L_mean":self.L_mean,
                "Li":self.Li,
                "Mc":self.Mc,
                "Mb":self.Mb,
                "age":self.age,}

        return (attributes)

    def import_attributes(self,attributes):
        self.DBH=attributes["DBH"]
        self.H=attributes["H"]
        self.CL=attributes["CL"]
        self.CD=attributes["CD"]
        self.CA=attributes["CA"]
        self.AGB=attributes["AGB"]
        self.LAI=attributes["LAI"]
        self.lmax=attributes["lmax"]
        self.lmin=attributes["lmin"]
        self.L_mean=attributes["L_mean"]
        self.Li=attributes["Li"]
        self.Mc=attributes["Mc"]
        self.Mb=attributes["Mb"]
        self.age=attributes["age"]
        self.Md=0
        #self.age+=1





    def H_from_DBH(self):
        """Calculate tree height calculated according to eq.1 in SI-Fisher et al.(2015).

        Returns:
            (float): H.
        """
        return self.h0*self.DBH**self.h1

    def CL_from_H(self):
        """Calculate tree crown_length (CL) calculated according to eq.2 in SI-Fisher et al.(2015)

        Returns:
            (float): CL.
        """
        return self.cl0*self.H

    def CD_from_DBH(self):
        """Calculate tree crown diameter (CD) calculated according to eq.3 in SI-Fisher et al.(2015).

        Returns:
            (float): CD.
        """
        return (self.cd0*self.DBH**self.cd1)-self.cd2

    def CA_from_CD(self):
        """Calculate tree crown area (CA) calculated according to eq.4 in SI-Fisher et al.(2015)

        Returns:
            (float): CA.
        """
        return (np.pi/4)*(self.CD**2)

    def calculate_volume(self):
        """ Calculate the volume of the stem based on DBH.

        Returns:
            (float): volume.
        """

        return self.H*(np.pi*((self.DBH/100)/2)**2)


    def AGB_from_DBH(self,D):
        """ Calculate tree aboveground biomass (AGB) calculated according to eq.5 in SI-Fisher et al.(2015)

        Returns:
            (float): AGB.
        """

        return (np.pi/4)*(D**2)*self.H*self.f*(self.rho/self.sigma)

    def calculate_f(self):
        """Calculate the form factor (f) calculated according to eq.6 in SI-Fisher et al.(2015)

        Returns:
            (float): f.

        """

        return self.f0*self.DBH**self.f1

    def LAI_from_DBH(self):
        """Calculate the leaf area index (LAI) calculated according to eq.7 in SI-Fisher et al.(2015)

        Returns:
            (float): LAI.
        """

        return self.l0*self.DBH**self.l1

    def Lind(self):
        """Calculate the incoming radiation on top of lmax layer the tree is reaching according to eq.36 in SI-Fisher et al.(2015).

        Returns:
            (float): Incoming radiation on top of this tree.
        """

        #crown_layers=(np.ceil(self.lmax)-np.floor(self.lmin))/self.world.topology.delta_h
        patch_LI=self.world.topology.LAIs[self.patch]
        top_layer=int(np.ceil(self.lmax))#/self.world.topology.delta_h)
        patch_LI=patch_LI[top_layer:]

        #self.patch_LI=sum(patch_LI)

        return self.world.topology.I0*np.exp(-self.k*sum(patch_LI))


#    def Lleaf(self,L):
#        """
#        Return the incoming radiation on top of leaves in layerL, according to eq.38 in SI-Fisher et al.(2015)
#        """

#        return (self.k/(1-self.m))*self.Li*np.exp(-self.k*L)

#    def Pleaf(self,L):
#        """
#        Return the gross photosynthetic rate on top of leaves in layer L, according to eq.37 in SI-Fisher et al.(2015)
#        """
#        lleaf=self.Lleaf(L)

#        return (self.alpha*lleaf*self.pmax)/(self.alpha*lleaf+pmax)

    def Pind(self,Li):
        """ Calculate the interim gross photosynthetic rate of one tree per year, according to eq.40 in SI-Fisher et al.(2015)

        Returns:
            (float) Gross photosynthetic rate.
        """

        pmk=self.pmax/self.k
        aki=self.alpha*self.k*Li
        pm=self.pmax*(1-self.m)

        pind=pmk*np.log((aki+pm)/(aki*np.exp(-self.k*self.LAI)+pm))

        #Tons of organic dry meter per year (todm.y-1)
        return pind*self.CA*60*60*self.lday*self.phi_act*np.power(10.0,-12)*2.27*44


    def calculate_rm(self):
            """Calculate the maintenance respiration rate(rm) according to eq.46 in SI-Fisher et al.(2015).

            Returns:
                (float): rm.
            """

            max_Li=self.world.topology.I0

            rm=(1/self.AGB)*(self.Pind(max_Li)-(max(0,(self.AGB_from_DBH(D=(self.DBH/100)+self.growth(D=self.DBH/100.))-self.AGB))/(1-self.rg)))

            return rm

    def growth(self,D):
        """Calculate the yearly diameter growth according to eq.48 in SI-Fisher et al.(2015). Assumes full availabilty of resources

        Returns:
            (float):diameter growth.
        """

        return (self.alpha0*D*(1-(D/self.Dmax)))*np.exp(-self.alpha1*D)

    def calculate_alpha0(self):
        """Calculate the alpha0 growth paramenter, according equations described in sctions F-5 (page 24) of SI-Fisher et al.(2015).

        Returns:
            (float): alpha 0.

        """


        exp=np.exp((self.Dmax-2*(self.DdeltaDmax*self.Dmax))/(self.Dmax-(self.DdeltaDmax*self.Dmax)))
        alpha0=exp*self.Dmax*self.deltaDmax
        alpha0=alpha0/((self.Dmax-(self.DdeltaDmax*self.Dmax))*(self.DdeltaDmax*self.Dmax))

        return alpha0

    def calculate_alpha1(self):
        """Calculate alpha1 growth paramenter, according equations described in sctions F-5 (page 24) of SI-Fisher et al.(2015).

        Returns:
            (float): alpha1.
        """

        alpha1=self.Dmax-2*(self.DdeltaDmax*self.Dmax)
        alpha1=alpha1/(self.Dmax*(self.DdeltaDmax*self.Dmax)-(self.DdeltaDmax*self.Dmax)**2)

        return alpha1

    def Biomass_gain(self):
        """ Calculate Aboveground Biomass (AGB) gain according to eq.43 in SI-Fisher et al.(2015).

        Returns:
            (float): AGB.
        """

        self.GPP=self.Pind(self.Li)
        self.Rm=self.calculate_rm()*self.AGB
        #self.Rm=self.r0*self.AGB+self.r1*self.AGB**2+self.r2*self.AGB**3
        gain=(1-self.rg)*(self.GPP-self.Rm)
        self.AGB+=gain
        return gain

    def DBH_from_Biomass(self,B):
        """ Calculate diameter at breast height from above ground biomass.

        Returns:
            (float): the DBH calculated from AGB.
        """


        p=(np.pi/4)*self.H*self.f*(self.rho/self.sigma)
        #if p>0 and B>0:
        dbh=np.sqrt(B/p)
        if isnan(dbh) or dbh<=0: return 0 #if Biomass is negative, new DBH will be 0
        return dbh

    def increase_DBH(self):
        """Increase tree DBH according to gain in biomass.

        Note: If DBH decreases to zero, remove tree.

        Returns:
            None.
        """
        gain=self.Biomass_gain()
        if round(self.AGB,3) <=0 or gain<0:
            self.remove_me()
        else:
            #new_B=self.AGB+self.Biomass_gain()
            self.DBH=self.DBH_from_Biomass(self.AGB)*100
            if round(self.DBH,3)<=2: self.remove_me()
        #if self.age>150: self.remove_me()

        #new_B=self.AGB+self.Biomass_gain()
        #self.DBH=self.DBH_from_Biomass(new_B)*100
        #if self.DBH<=0: self.remove_me()

    def calculate_Rc(self):
        patch_CCA=self.world.topology.CCAs[self.patch]
        max_CCA=max(patch_CCA[int(np.ceil(self.lmin)):int(np.ceil(self.lmax))],default=self.CA)

        if max_CCA==0: return 1.0
        rc=1/max_CCA
        if rc>=1.0:
            return 0
        else:
            return (1-rc)

    def calculate_Ms(self):
        """Calculate the size-dependent mortality for small trees according to model description in Ruger et al. 2007. (SI, page 4).

        Returns:
            (float): Ms.
        """
        Ms=0


        if self.DBH <self.Dmort:
            Ms=self.m_max-self.m_max*(self.DBH/self.Dmort)
        #elif self.DBH>=self.max_dbh:
        #    Ms=self.m_max

        return Ms


    def stochastic_M(self):
        """Stochastic mortality.

        Returns:
            (bools): True if tree should be removed.
        """
        #if self.DBH>=self.max_dbh:
        #    return True #np.random.random()<= self.m_max
        w_Mb=0
        w_Mc=1
        w_Md=1


        M=((w_Mc*self.Mc+w_Mb*self.Mb+w_Md*self.Md)/(w_Mc+w_Mb+w_Md))#
        return np.random.random()<=M#*self.world.mortality_factor

    def stochastic_Cm(self):
        """Stochastic crownding mortality.

        Returns:
            (bool): True if tree dies from crowding.
        """
        return np.random.random()<=min(self.Mc,1.0)

    def stochastic_Bm(self):
        """Stochastic Background mortality.

        Returns:
            (bool): True if tree dies from background mortality.
        """
        return np.random.random()<=min(self.Mb,1.0)

    def stochastic_Dm(self):
        """Stochastic damage mortality.

        Returns:
            (bool): True if tree dies from damage.
        """
        return np.random.random()<=min(self.Md,1.0)

    def log_me(self):
        """Cut down this tree.

            Randomly determines where tree is going to fall and inflict damage in the the trees hit by the crown. Only trees with DBH< 50 cm are damaged.

            The carbon of this trees is not added to any Stock.
        """


        subclass="Tree_"+self.Ftype
        globals()[subclass].Indices.remove(self.id)
        self.world.topology.trees_per_patch[self.patch].remove(self.id)
        #self.world.Smort+=self.AGB*0.44

        topology=self.world.topology

        dir=np.random.randint(1,360)
        x=self.position[0]+self.H*np.sin(np.deg2rad(2*np.pi*(dir/360)))
        y=self.position[1]+self.H*np.cos(np.deg2rad(2*np.pi*(dir/360)))
        md=min(self.CA/topology.patch_area,1)

        target_patch=np.floor((x,y))
        target_patch=(int(target_patch[0]),int(target_patch[1]))
        target_patch=topology.in_bounds(target_patch)



        neighborhood=topology.hood(cell=target_patch,remove_center=False)
        trees_in_neighborhood=topology.trees_in_patches(neighborhood)
        trees_in_circle=[id for id,p in trees_in_neighborhood if topology.point_within_circle(x,y,(self.CD/2),p[0],p[1],topology.scaled_distance,scale=topology.patch_area**0.5)]

        for t in trees_in_circle:
            if Tree.Instances.get(t).DBH <=50:
                 Tree.Instances.get(t).Md=md



        del Tree.Instances[self.id]


    def remove_me(self):
        """
        Kill this tree.
        """
        subclass="Tree_"+self.Ftype
        globals()[subclass].Indices.remove(self.id)
        if self.DBH>self.Dfall and np.random.random()<=self.pfall:
            self.fall()
        self.world.topology.trees_per_patch[self.patch].remove(self.id)
        self.world.Smort+=max(0,self.AGB)*0.44
        Tree.DeadTrees[self.Ftype]+=1
        del Tree.Instances[self.id]

    def fall(self):
        """
        Inflict damage on other trees affected by the falling one.
        """
        topology=self.world.topology

        dir=np.random.randint(1,360)
        x=self.position[0]+self.H*np.sin(np.deg2rad(2*np.pi*(dir/360)))
        y=self.position[1]+self.H*np.cos(np.deg2rad(2*np.pi*(dir/360)))
        md=min(self.CA/topology.patch_area,1)

        #Damage caused on trees near the falling trees, due to lianas

        neighborhood=topology.hood(cell=self.patch,remove_center=False)
        trees_in_neighborhood=topology.trees_in_patches(neighborhood)
        trees_in_circle=[id for id,p in trees_in_neighborhood if topology.point_within_circle(self.position[0],self.position[1],(self.CD/2),p[0],p[1],topology.scaled_distance,scale=topology.patch_area**0.5)]

        for t in trees_in_circle:
            if Tree.Instances.get(t).H <=(self.H+1):
                 Tree.Instances.get(t).Md=md*(topology.surface[topology.dim_names["LianaCoverage"]][self.patch])



        target_patch=np.floor((x,y))
        target_patch=(int(target_patch[0]),int(target_patch[1]))
        target_patch=topology.in_bounds(target_patch)



        neighborhood=topology.hood(cell=target_patch,remove_center=False)
        trees_in_neighborhood=topology.trees_in_patches(neighborhood)
        trees_in_circle=[id for id,p in trees_in_neighborhood if topology.point_within_circle(x,y,(self.CD/2),p[0],p[1],topology.scaled_distance,scale=topology.patch_area**0.5)]

        for t in trees_in_circle:
            if Tree.Instances.get(t).H <=(self.H+1):
                 Tree.Instances.get(t).Md=md




###############################################################################
#                                 Tree Types                                  # ###############################################################################


class Tree_FT1(Tree):
    #Instances={}
    Indices=[]
    wood_density=0.35
    max_agb=200
    max_dbh=145.196
    h0=3.28
    h1=0.57
    cl0=0.80
    cd0=0.60
    cd1=0.68
    cd2=0.0
    rho=0.55
    sigma=0.70
    f0=0.77
    f1=-0.18
    l0=2.0
    l1=0.10
    Iseed=0.03
    Nseed=30
    m=0.5
    rg=0.25
    Mb=0.015
    DdeltaDmax=0.33
    deltaDmax=0.012
    Ftype='FT1'
    pmax=2.0
    alpha=0.36



    #def __repr__(self):
    #    return "Functional Type 1-id: {0}".format(self.id)


    def __init__(self,position,world,dbh,age=0,id=None):
        super().__init__(position=position,world=world,dbh=dbh,id=id,wood_density=Tree_FT1.wood_density,max_agb=Tree_FT1.max_agb,max_dbh=Tree_FT1.max_dbh,h0=Tree_FT1.h0,h1=Tree_FT1.h1,cl0=Tree_FT1.cl0,cd0=Tree_FT1.cd0,cd1=Tree_FT1.cd1,cd2=Tree_FT1.cd2,rho=Tree_FT1.rho,sigma=Tree_FT1.sigma,f0=Tree_FT1.f0,f1=Tree_FT1.f1,l0=Tree_FT1.l0,l1=Tree_FT1.l1,m=Tree_FT1.m,rg=Tree_FT1.rg,Mb=world.mortality_factor*Tree_FT1.Mb,DdeltaDmax=Tree_FT1.DdeltaDmax,deltaDmax=Tree_FT1.deltaDmax,Ftype=Tree_FT1.Ftype,pmax=Tree_FT1.pmax,alpha=Tree_FT1.alpha,age=age)



        if id == None:
            #self.id=Tree.ID
            Tree_FT1.Indices.append(self.id)
            Tree.Instances[self.id]=self
            Tree.IncrementID()
        else:
            self.id=id
            Tree_FT1.Indices.append(self.id)
            Tree.Instances[self.id]=self





class Tree_FT2(Tree):
    Indices=[]
    wood_density=0.35
    max_agb=200
    max_dbh=120
    h0=4.64
    h1=0.41
    cl0=0.80
    cd0=0.60
    cd1=0.68
    cd2=0.0
    rho=0.54
    sigma=0.70
    f0=0.77
    f1=-0.18
    l0=2.0
    l1=0.10
    Iseed=0.01
    Nseed=156
    m=0.5
    rg=0.25
    Mb=0.03
    pmax=3.1
    alpha=0.28
    DdeltaDmax=0.34
    deltaDmax=0.012
    Ftype='FT2'

    #Instances={}



    def __init__(self,position,world,dbh,age=0,id=None):
        super().__init__(position=position,world=world,dbh=dbh,id=id,wood_density=Tree_FT2.wood_density,max_agb=Tree_FT2.max_agb,max_dbh=Tree_FT2.max_dbh,h0=Tree_FT2.h0,h1=Tree_FT2.h1,cl0=Tree_FT2.cl0,cd0=Tree_FT2.cd0,cd1=Tree_FT2.cd1,cd2=Tree_FT2.cd2,rho=Tree_FT2.rho,sigma=Tree_FT2.sigma,f0=Tree_FT2.f0,f1=Tree_FT2.f1,l0=Tree_FT2.l0,l1=Tree_FT2.l1,m=Tree_FT2.m,rg=Tree_FT2.rg,Mb=world.mortality_factor*Tree_FT2.Mb,DdeltaDmax=Tree_FT2.DdeltaDmax,deltaDmax=Tree_FT2.deltaDmax,Ftype=Tree_FT2.Ftype,pmax=Tree_FT2.pmax,alpha=Tree_FT2.alpha,age=age)


        if id == None:
            #self.id=Tree.ID
            Tree_FT2.Indices.append(self.id)
            Tree.Instances[self.id]=self
            Tree.IncrementID()
        else:
            self.id=id
            Tree_FT2.Indices.append(self.id)
            Tree.Instances[self.id]=self




class Tree_FT3(Tree):
    Indices=[]
    wood_density=0.35
    max_agb=200
    max_dbh=80
    h0=4.82
    h1=0.44
    cl0=0.30
    cd0=0.60
    cd1=0.68
    cd2=0.0
    rho=0.41
    sigma=0.70
    f0=0.77
    f1=-0.18
    l0=2.0
    l1=0.10
    Iseed=0.05
    Nseed=21
    m=0.5
    rg=0.25
    Mb=0.029
    pmax=6.8
    alpha=0.23
    DdeltaDmax=0.23
    deltaDmax=0.019
    Ftype='FT3'



    def __init__(self,position,world,dbh,age=0,id=None):
        super().__init__(position=position,world=world,dbh=dbh,id=id,wood_density=Tree_FT3.wood_density,max_agb=Tree_FT3.max_agb,max_dbh=Tree_FT3.max_dbh,h0=Tree_FT3.h0,h1=Tree_FT3.h1,cl0=Tree_FT3.cl0,cd0=Tree_FT3.cd0,cd1=Tree_FT3.cd1,cd2=Tree_FT3.cd2,rho=Tree_FT3.rho,sigma=Tree_FT3.sigma,f0=Tree_FT3.f0,f1=Tree_FT3.f1,l0=Tree_FT3.l0,l1=Tree_FT3.l1,m=Tree_FT3.m,rg=Tree_FT3.rg,Mb=world.mortality_factor*Tree_FT3.Mb,DdeltaDmax=Tree_FT3.DdeltaDmax,deltaDmax=Tree_FT3.deltaDmax,Ftype=Tree_FT3.Ftype,pmax=Tree_FT3.pmax,alpha=Tree_FT3.alpha,age=age
        )

        if id == None:
            #self.id=Tree.ID
            Tree_FT3.Indices.append(self.id)
            Tree.Instances[self.id]=self
            Tree.IncrementID()
        else:
            self.id=id
            Tree_FT3.Indices.append(self.id)
            Tree.Instances[self.id]=self


class Tree_FT4(Tree):
    Indices=[]
    wood_density=0.35
    max_agb=200
    max_dbh=80
    h0=4.27
    h1=0.43
    cl0=0.30
    cd0=0.60
    cd1=0.68
    cd2=0.0
    rho=0.40
    sigma=0.70
    f0=0.77
    f1=-0.18
    l0=2.0
    l1=0.10
    Iseed=0.02
    Nseed=300
    m=0.5
    rg=0.25
    Mb=0.04
    pmax=11
    alpha=0.20
    DdeltaDmax=0.60
    deltaDmax=0.029
    Ftype='FT4'



    def __init__(self,position,world,dbh,age=0,id=None):
        super().__init__(position=position,world=world,dbh=dbh,id=id,wood_density=Tree_FT4.wood_density,max_agb=Tree_FT4.max_agb,max_dbh=Tree_FT4.max_dbh,h0=Tree_FT4.h0,h1=Tree_FT4.h1,cl0=Tree_FT4.cl0,cd0=Tree_FT4.cd0,cd1=Tree_FT4.cd1,cd2=Tree_FT4.cd2,rho=Tree_FT4.rho,sigma=Tree_FT4.sigma,f0=Tree_FT4.f0,f1=Tree_FT4.f1,l0=Tree_FT4.l0,l1=Tree_FT4.l1,m=Tree_FT4.m,rg=Tree_FT4.rg,Mb=world.mortality_factor*Tree_FT4.Mb,DdeltaDmax=Tree_FT4.DdeltaDmax,deltaDmax=Tree_FT4.deltaDmax,Ftype=Tree_FT4.Ftype,pmax=Tree_FT4.pmax,alpha=Tree_FT4.alpha,age=age)


        if id == None:
            #self.id=Tree.ID
            Tree_FT4.Indices.append(self.id)
            Tree.Instances[self.id]=self
            Tree.IncrementID()
        else:
            self.id=id
            Tree_FT4.Indices.append(self.id)
            Tree.Instances[self.id]=self

class Tree_FT5(Tree):
    Indices=[]
    wood_density=0.35
    max_agb=200
    max_dbh=47
    h0=4.35
    h1=0.34
    cl0=0.30
    cd0=0.60
    cd1=0.68
    cd2=0.0
    rho=0.52
    sigma=0.70
    f0=0.77
    f1=-0.18
    l0=2.0
    l1=0.10
    Iseed=0.03
    Nseed=2
    m=0.5
    rg=0.25
    Mb=0.021
    pmax=7
    alpha=0.30
    DdeltaDmax=0.33
    deltaDmax=0.011
    Ftype='FT5'



    def __init__(self,position,world,dbh,age=0,id=None):
        super().__init__(position=position,world=world,dbh=dbh,id=id,wood_density=Tree_FT5.wood_density,max_agb=Tree_FT5.max_agb,max_dbh=Tree_FT5.max_dbh,h0=Tree_FT5.h0,h1=Tree_FT5.h1,cl0=Tree_FT5.cl0,cd0=Tree_FT5.cd0,cd1=Tree_FT5.cd1,cd2=Tree_FT5.cd2,rho=Tree_FT5.rho,sigma=Tree_FT5.sigma,f0=Tree_FT5.f0,f1=Tree_FT5.f1,l0=Tree_FT5.l0,l1=Tree_FT5.l1,m=Tree_FT5.m,rg=Tree_FT5.rg,Mb=world.mortality_factor*Tree_FT5.Mb,DdeltaDmax=Tree_FT5.DdeltaDmax,deltaDmax=Tree_FT5.deltaDmax,Ftype=Tree_FT5.Ftype,pmax=Tree_FT5.pmax,alpha=Tree_FT5.alpha,age=age)


        if id == None:
            self.id=Tree.ID
            Tree_FT5.Indices.append(self.id)
            Tree.Instances[self.id]=self
            Tree.IncrementID()
        else:
            self.id=id
            Tree_FT5.Indices.append(self.id)
            Tree.Instances[self.id]=self


class Tree_FT6(Tree):
    Indices=[]
    wood_density=0.35
    max_agb=200
    max_dbh=16.123
    h0=3.0
    h1=0.60
    cl0=0.30
    cd0=0.60
    cd1=0.68
    cd2=0.0
    rho=0.47
    sigma=0.70
    f0=0.77
    f1=-0.18
    l0=2.0
    l1=0.10
    Iseed=0.02
    Nseed=200
    m=0.5
    rg=0.25
    Mb=0.045
    pmax=12
    alpha=0.20
    DdeltaDmax=0.60
    deltaDmax=0.029
    Ftype='FT6'


    def __init__(self,position,world,dbh,age=0,id=None):
        super().__init__(position=position,world=world,dbh=dbh,id=id,wood_density=Tree_FT6.wood_density,max_agb=Tree_FT6.max_agb,max_dbh=Tree_FT6.max_dbh,h0=Tree_FT6.h0,h1=Tree_FT6.h1,cl0=Tree_FT6.cl0,cd0=Tree_FT6.cd0,cd1=Tree_FT6.cd1,cd2=Tree_FT6.cd2,rho=Tree_FT6.rho,sigma=Tree_FT6.sigma,f0=Tree_FT6.f0,f1=Tree_FT6.f1,l0=Tree_FT6.l0,l1=Tree_FT6.l1,m=Tree_FT6.m,rg=Tree_FT6.rg,Mb=world.mortality_factor*Tree_FT6.Mb,DdeltaDmax=Tree_FT6.DdeltaDmax,deltaDmax=Tree_FT6.deltaDmax,Ftype=Tree_FT6.Ftype,pmax=Tree_FT6.pmax,alpha=Tree_FT6.alpha,age=age)



        if id == None:
            self.id=Tree.ID
            Tree_FT6.Indices.append(self.id)
            Tree.Instances[self.id]=self
            Tree.IncrementID()
        else:
            self.id=id
            Tree_FT6.Indices.append(self.id)
            Tree.Instaces[self.id]=self
