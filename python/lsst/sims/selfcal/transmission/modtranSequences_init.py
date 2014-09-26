#####
#  Class used to generate sequences of atmospheric parameters for input to MODTRAN.
###

import numpy

class AtmosphereSequences():
    """
    This class is used to generate sequences of atmospheric parameters, matching long-term weather trends,
    and accounting for the airmass/pointing location of each observation. These atmospheric parameters are the
    basic inputs for MODTRAN, so can be used to generate an atmospheric transmission curve. The class generates
    these parameters for a series of opsim pointings.
    """
    
    def __init__(self):
        """Instantiate an AtmosphereSequences object."""
        self.modtran_visits = []
        # individual modtran visit = dictionary with necessary keys
        self.modtran_keys = ('O3', 'H2O', 'IHAZE', 'ISEAS', 'IVULC', 'ICLD', 'IVSA', 'VIS',
                             'ZAER11', 'ZAER12', 'SCALE1', 'ZAER21', 'ZAER22', 'SCALE2')
        """
        # so, for example ..
        self.modtran_visits[i] = {}
        for key in self.modtran_keys:
            self.modtran_visits[i][key] = value
        """
        return


    def generate_parameters(self, opsim_visits):
        """Generate the atmospheric parameters over time.

        The input parameter, 'opsim_visits', is a dictionary of numpy arrays containing opsim visit information
        as opsim_visits['expmjd'] (a numpy array of expmjds), 'ra', 'dec', and 'obshistid' (at a minimum).
        These arrays are already sorted order of increasing expmjd. 
        Returns a list of dictionaries containing the modtran information for each opsim visit (in the same order)."""
        # I'm not sure if the parameter history is first generated for the atmosphere at zenith, and then adjusted
        # for the location of the actual field or not, and if there are different timescales for this generation.
        # Obviously, the methods needed here will change a little bit as a result.
        #   (perhaps there will be a 'generate_preliminary sequence', then a generate exact parameters method..
        #    ... and maybe even individual methods for each component of the atmosphere ? )
        # Eventually, there should be a LIST of DICTIONARIES which holds all the atmospheric parameters which
        # need to be input into modtran.
        # Each dictionary holds the 12 parameters ('O3', 'H2O', 'ihaze', ...) and then there will be a list of these
        # dictionaries to hold the series of information over time. 
        #  In the end we will write these parameters into a database. Here, all we have to do is
        # generate the parameters. Please add other methods as needed .... (is there is one per atmosphere component?)
        # MUST make sure not to shuffle the obshistid / expmjd list - return list of MODTRAN parameters in the same order.
        
        return self.modtran_visits


