from py_ibm import *
from trees_ibm.tree_agent import Tree
#from plot_trees import plot_tree_pos
import pdb
from tables import *
import json
import numpy as np


class Tree_Grid(Rectangular_Grid):
    """ Represent space as a grid.


    Args:
        x_max : int
            Number of horizontal cells
        y_max : int
            Number of vertical cells
        patch_area : float
            Area of each patch (cell) in squared meters
        I0 : float
            Incoming irradiance on top of canopy [umol photon / m^2s]
        k : float
            Light extinction coefficient.
        delta_h :float
            Height of each vertical layer.
        phi_act : float
            Photosynthetic active period.
        lday : int
            Mean day length during the vegetation period phi_act.
        ndim : int
            Number of dimensions. One per environmental variable.
        dim_names : list
            A list of 'ndim' strings representing the names of dimensions.

    Attributes:
        LAIs : dict
            Dictionary with patches as keys and a list of averaged Leaf Area Indices (one for each vertical layer) as values.
        CCAs : dict
            Dictionary with patches as keys a a list of Cummulative Crown Areas (one for each vertical layer) as values.
        trees_per_patch : list
            Dictionary with patches ('(x,y)') as keys and a list of tree ids as values
        total_area :float
            Area of all simulated patches in ha
    """

    def __init__(self, x_max, y_max, patch_area, I0, k, delta_h, phi_act, lday, ndim=5, dim_names=["GroundLight", "LianaCoverage", "Temperature", "WindSpeed", "Humidity"]):
        super().__init__(x_max, y_max, ndim=ndim, dim_names=dim_names)
        self.patch_area = patch_area
        self.delta_h = delta_h
        self.I0 = I0
        self.k = k
        self.lday = lday
        self.phi_act = phi_act
        self.LAIs = {}
        self.CCAs = {}
        self.total_area = self.x_max * self.y_max * self.patch_area / 10000
        # s={}
        self.trees_per_patch = {(x, y): [] for x in range(
            0, self.x_max) for y in range(0, self.y_max)}
        # self.FTs=["FT1","FT2","FT3","FT4","FT5","FT6"]
        self.FTs = [k for k in Tree.DERIVED_TREES.keys()]
        self.seedbank = {k: [] for k in self.FTs}

    def exeeding_CCA(self):
        """ Finds which patches have the cummulative crown area bigger than the patch area in at least one layer.

            Returns: list
                Alist of patchs ([(x,y),(x,y),...,(x,y)])
        """
        patches = [patch for patch, layers in self.CCAs.items() if any(
            layer > self.patch_area for layer in layers)]

        return patches

    def update(self):
        """ Execute the methods that calculate values for each dimension of the grid.

        Returns:
            None.

        """

        self.LAIs = self.patch_based_LAI(Tree.TreesPerPatch())
        self.CCAs = self.patch_based_CCA(Tree.TreesPerPatch())
        self.ground_level_light()
        self.surface[self.dim_names["LianaCoverage"]
                     ] = self.calculate_liana_coverage()

    def seed_establishment(self, seedbank):
        """Check if the seeds in the seedbank have the light and space conditions to germinate.

            Note: All seeds in the given seedbank are assumed to have the same conditions. Multiple calls to this method passing different seedbanks allow the representation of different requirements.

            All established individuals are initiated with a predefined DBH (the same for all individuals), but initial height is calculated based on functional relationships and only if there is enough space in that layer a new individual will be stablished.

        Args:

            seedbank : dict
                Dictionary with Functional Type (FT) and list of seed positions as values.
            lmin : float
                Minimum light intensity at ground level for seed to germinate.
            h0 : float
                Allometric paramenter used to calculate initial height.
            h1 : float
                Allometric paramenter used to calculate initial height.

        Returns: dict
            est_seeds: a dictionary with patches as keys ('(x,y)') and number of established seeds as values.
        """
        est_seeds = {}

        for ft in seedbank.keys():
            h0 = Tree.DERIVED_TREES[ft].h0
            h1 = Tree.DERIVED_TREES[ft].h1
            min_light = Tree.DERIVED_TREES[ft].Iseed
            seeds_to_keep = []
            for i in range(len(seedbank[ft])):
                #print(i, len(seedbank[ft]), ft)
                seed = seedbank[ft][i]
                patch = (int(np.floor(seed[0])), int(np.floor(seed[1])))
                if np.random.rand() > 0.9995:
                    seeds_to_keep.append(i)
                    if not self.check_light(min_light, patch):
                        seeds_to_keep.remove(i)
                        # break
                    elif not self.check_CCA(h0, h1, patch):
                        seeds_to_keep.remove(i)

            est_seeds[ft] = [seedbank[ft][s] for s in seeds_to_keep]
            # est_seeds[ft]=seeds_to_keep

        return est_seeds

    def check_light(self, min_light, patch):
        """ Check if the light intensity at ground level meets the needs for seed germination.

            Note: Ligth intensity is homogeneized at the patch level.

        Args:
            min_light : float
                Minimum light intensity at ground for seeds to germinate.
            patch : tuple
                (x,y), patch  that will have the light conditions assessed.

        Returns: bool.
            True if light level is above minimum.
        """
        return self.surface[self.dim_names["GroundLight"]][patch] > min_light

    def check_CCA(self, h0, h1, patch):
        """ Check if the there is enough space forthe new plants to be established

            Note: Cummulative Crown Area is homogeneized at the patch level

        Args:
            h0 : float
                Allometric paramenter used to calculate initial height
            h1 : float
                Allometric paramenter used to calculate initial height
            patch : tuple
                (x,y), patch  that will have the light conditions assessed.

        Returns: bool.
            True if there is enough space.
        """
        Hmin = h0 * 2**h1
        l = int(np.floor(Hmin / self.delta_h))
        if self.CCAs.get(patch) == None:
            return True
        else:
            return self.CCAs[patch][l] < 1.0

    def ground_level_light(self):
        """Calculate the amount of light that reaches the ground.

            The value per patch is stored in the "GroundLight" dimension.
            Values are homogeneized at the patch level.

            Returns:
                None.
        """
        self.surface[self.dim_names["GroundLight"]] = self.I0
        for patch in self.LAIs.keys():
            self.surface[self.dim_names["GroundLight"]][patch] = self.I0 * \
                np.exp(-self.k * sum(self.LAIs[patch]))

    def patch_based_LAI(self, trees_per_patch):
        """
        Calculate the Leaf Area Index for each patch in the grid.


        Args:
            trees_per_patch :dict
                Dictionary in which keys are patch coordinates (x,y) and values are a list of the tree ids in that patch.

        Returns: dict
            L_per_layer: dictionary in which the keys are patch coordinates('(x,y)')-the same as the ones in trees_per_patch- and values are a list of the LAI per vertical layer
        """

        L_per_layer = {}
        nlayers = int(np.ceil(np.ceil(Tree.MaxHeight()) / self.delta_h))
        for k in trees_per_patch.keys():
            L_per_layer[k] = []
            lower_limit = int(
                min([Tree.Instances.get(t).lmin for t in trees_per_patch[k]]))
            upper_limit = int(
                max([Tree.Instances.get(t).lmax for t in trees_per_patch[k]]))

            for i in range(nlayers + 1):

                L_per_layer[k].append((1 / self.patch_area) * sum([Tree.Instances.get(
                    tree_id).L_mean for tree_id in trees_per_patch[k] if lower_limit <= i and upper_limit >= i]))

                # upper and lower limits
                #Tree.Instances.get(tree_id).lmin<=i*self.delta_h and Tree.Instances.get(tree_id).lmax>=i*self.delta_h
        return L_per_layer

    def patch_based_CCA(self, trees_per_patch):
        """Calculate the Leaf Area Index for each patch in the grid.

        Args:
            trees_per_patch : dict
            Dictionary in which keys are patches ('(x,y)')' and values are lists of the trees ids in that patch.

        Returns: dict
            L_per_layer: dictionary in which the keys are the same as the ones in trees_per_patch and values are a list of the CCA per vertical layer

            """
        L_per_layer = {}
        nlayers = int(np.ceil(np.ceil(Tree.MaxHeight()) / self.delta_h))
        for k in trees_per_patch.keys():
            L_per_layer[k] = []
            lower_limit = int(
                min([Tree.Instances.get(t).lmin for t in trees_per_patch[k]]))
            upper_limit = int(
                max([Tree.Instances.get(t).lmax for t in trees_per_patch[k]]))

            for i in range(nlayers + 1):

                L_per_layer[k].append((1 / self.patch_area) * sum([Tree.Instances.get(
                    tree_id).CA for tree_id in trees_per_patch[k] if lower_limit <= i and upper_limit >= i]))
        return L_per_layer

    def trees_in_patches(self, patches):
        """ Make a list of tree ids and their respective positions.

        Args:
            patch : list
                A list of tuples representing the patch coordinates '(x,y)'.

        Returns: list
            id_pos: a list tuples with tree ids and positions '(id,(x_pos,y_pos))'

        """

        # one list of ids for each patch
        ids = [self.trees_per_patch[p] for p in patches]
        # one single list of tuples with (id, pos)
        id_pos = [(tree, Tree.Instances.get(tree).position)
                  for patch in ids for tree in patch]

        return id_pos

    def trees_per_age_groups(self, list_of_trees):
        """
        Return a dictionary in which keys are a age group and values are lists containing the id of the respective trees.

            Args:
                list_of_trees: list
                    A list with the id of trees.

            Returns: dict
                A dictionary age as keys and a list of tree ids as values.

        """

        age_groups = {}
        for t in list_of_trees:
            if age_groups.get(Tree.Instances[t].age) == None:
                age_groups[Tree.Instances[t].age] = [Tree.Instances[t].id]
            else:
                age_groups[Tree.Instances[t].age].append(Tree.Instances[t].id)

        return age_groups

    def trees_per_type(self, list_of_trees):
        """Organize trees by functional type.

        Args:
            list_of_trees :list
                A list of tree ids.

        Returns: dict
            type_groups: a dictionary in which keys are a Functional
            Types and values are lists containing the id of the respective trees.
        """

        type_groups = {}
        for t in list_of_trees:
            if type_groups.get(Tree.Instances[t].Ftype) == None:
                type_groups[Tree.Instances[t].Ftype] = [Tree.Instances[t].id]
            else:
                type_groups[Tree.Instances[t].Ftype].append(
                    Tree.Instances[t].id)

        return type_groups

    def cohorts(self, patches):
        """ Assemble cohorts based on tree type and age.

            Uses the trees_per_type() and trees_by_age_groups() methods to
            organize a dictionary of dictionaries containing the indices of
            trees in each age group within each type

        Args:
            Patches :list
                A list of tuples representing patch coordinatas
            [(x1,y1),(x2,y2)...]

        Returns: dict
            cohorsts: a dictionary with trees organized by type and age:
             {FunctionalType1:{AgeGroup1:[tree1,tree2,tree3...]}}
        """

        trees_in_patch = self.trees_in_patches(patches)
        list_of_trees = [i[0] for i in trees_in_patch]
        FT = self.trees_per_type(list_of_trees)
        cohorts = {}
        for group in FT.keys():
            cohorts[group] = self.trees_per_age_groups(FT[group])

        return cohorts

    def calculate_liana_coverage(self):
        """ Calculate the coverage (%) os lianas in on patch based on the previous coverage, temperature, humidity and wind speed

        Returns:
                None
        """
        return np.zeros((self.x_max, self.y_max))

    def calculate_liana_biomass(self):
        """
        Calculate the biomass of lianas in one patch based on the coverage.
        """
        pass
