import numpy as np

class SLR_C20(object):
    '''
    class SLR_C20
    - attributes:
        - info 
        - normalization
        - background_gravity
        - title
        - summary
        - institution
        - product_version
        - time_coverage_start
        - time_coverage_end
        - total_month
        - total_month_counts
        - solution_month
        - solution_counts
        - missing_month
        - missing_month_counts
        - missing_solution_flag
        - date_issued
        - mean_c20
        - c20
        - c20_demean
        - c20_std
    
    - methods:
        - deaverage
    '''
    def __init__(self,info,c20,c20_demean,c20_std):
        self.info = info
        self.normalization = info['normalization'] 
        self.background_gravity = info['background_gravity']
        self.title = info['title']
        self.summary = info['summary']
        self.institution = info['institution']
        self.product_version = info['product_version']
        self.time_coverage_start = info['time_coverage_start']
        self.time_coverage_end = info['time_coverage_end']
        self.total_month = info['total_month']
        self.total_month_counts = info['total_month_counts'] 
        self.solution_month = info['solution_month']
        self.solution_counts = info['solution_counts'] 
        self.missing_month = info['missing_month']
        self.missing_month_counts = info['missing_month_counts']
        self.missing_solution_flag = info['missing_solution_flag']
        self.date_issued = info['date_issued']
        self.mean_c20 = info['mean_c20']
        
        self.c20 = c20
        self.c20_demean = c20_demean
        self.c20_std = c20_std
        
    def __repr__(self):
    
        return 'title = {:s}\nnormalization = {:s}\ninstitution = {:s}\nproduct_version = {:s}\ntime_coverage_start = {:s}\ntime_coverage_end = {:s}\nsolution_counts = {:d}\ntotal_month_counts = {:d}\nmissing_month_counts = {:d}'.format\
        (self.title,self.normalization,self.institution,self.product_version,self.time_coverage_start,self.time_coverage_end,self.solution_counts,self.total_month_counts,self.missing_month_counts)
        
    def deaverage(self):
        
        '''
        Deaverage the SLR C20 data.
    
        Usage: 
        slr_c20_deaverage = slr_c20.deaverage()

        Outputs:
        slr_c20_deaverage -> Deaveraged SLR C20 solitions
        
        Examples:
        -----------
        >>> slr_c20 = read_SLR_C20_file()
        >>> slr_c20_deaverage = slr_c20.deaverage()
        >>> print(slr_c20_deaverage)
        [ 1.58283192e-10,  1.20785392e-10, -4.92720085e-11, -3.20084085e-11,
         -1.73326085e-11,  9.12899151e-12,  5.35746915e-11,  9.35680915e-11,
          8.24645915e-11,  1.23541892e-10,  1.03958091e-10,  9.27922915e-11,
          ...
          ...
         -1.44952909e-10, -1.42519708e-10, -9.05243085e-11, -7.39151085e-11,
         -6.21788085e-11, -2.70057085e-11, -4.02519085e-11,  1.70110915e-11,
          2.92669150e-12, -4.49227085e-11, -1.89951408e-10, -1.86842409e-10]
        '''
        info = self.info.copy()

        background_mean = np.average(self.c20)

        info['background_gravity'] = 'Average of monthly solutions'   
        info['mean_c20'] = background_mean
        c20_deaverage = self.c20 - background_mean
        c20_demean = c20_deaverage
        
        info['title'] = 'Deaveraged ' + info['title']
        
        return SLR_C20(info,c20_deaverage,c20_demean,self.c20_std)

    def __add__(self,other):

        info = self.info.copy()
        #info['title'] = 
        #info['summary'] =

        self_c20,self_c20_demean,self_c20_std = self.c20,self.c20_demean,self.c20_std 
        other_c20,other_c20_demean,other_c20_std = other.c20,other.c20_demean,other.c20_std

        if 'rate' in self.title and 'rate' in other.title:
            add_c20 = self_c20 + other_c20
            add_c20_demean = self_c20_demean + other_c20_demean
            add_c20_std = np.sqrt(self_c20_std**2 + other_c20_std**2)

        elif ('rate' in self.title and 'rate' not in other.title) or ('rate' not in self.title and 'rate' in other.title):
            raise Exception('addition can not be completed between rate object and series object') 

        else: 
            self_existing_solution_flag = ~self.missing_solution_flag
            other_existing_solution_flag = ~other.missing_solution_flag

            if (self_existing_solution_flag != other_existing_solution_flag).any():
                existing_solution_flag = self_existing_solution_flag & other_existing_solution_flag

                n = len(existing_solution_flag)

                assumed_self_c20 = np.zeros(n)
                assumed_other_c20 = assumed_self_c20.copy()

                assumed_self_c20[self_existing_solution_flag] = self.c20
                assumed_other_c20[other_existing_solution_flag] = other.c20

                self_c20 = assumed_self_c20[existing_solution_flag] 
                other_c20 = assumed_other_c20[existing_solution_flag] 

                assumed_self_c20_demean = np.zeros(n)
                assumed_other_c20_demean = assumed_self_c20_demean.copy()

                assumed_self_c20_demean[self_existing_solution_flag] = self.c20_demean
                assumed_other_c20_demean[other_existing_solution_flag] = other.c20_demean

                self_c20_demean = assumed_self_c20_demean[existing_solution_flag] 
                other_c20_demean = assumed_other_C20_demean[existing_solution_flag] 

                assumed_self_c20_std = np.zeros(n)
                assumed_other_c20_std = assumed_self_c20_std.copy()

                assumed_self_c20_std[self_existing_solution_flag] = self.c20_std
                assumed_other_c20_std[other_existing_solution_flag] = other.c20_std

                self_c20_std = assumed_self_c20_std[existing_solution_flag] 
                other_c20_std = assumed_other_c20_std[existing_solution_flag] 

                solution_month = self.total_month[existing_solution_flag]
                solution_counts = len(solution_month)
                missing_solution_flag = ~existing_solution_flag
                missing_month = self.total_month[missing_solution_flag]
                missing_month_counts = len(missing_month)

                info['solution_month'] = solution_month
                info['solution_counts'] = solution_counts
                info['missing_month'] = missing_month
                info['missing_month_counts'] = missing_month_counts
                info['missing_solution_flag'] = missing_solution_flag
        
            add_c20 = self_c20 + other_c20
            add_c20_demean = self_c20_demean + other_c20_demean
            add_c20_std = np.sqrt(self_c20_std**2 + other_c20_std**2)

        return SLR_C20(info,add_c20,add_c20_demean,add_c20_std)
    
    def __sub__(self,other):

        info = self.info.copy()
        #info['title'] = 
        #info['summary'] =

        self_c20,self_c20_demean,self_c20_std = self.c20,self.c20_demean,self.c20_std
        other_c20,other_c20_demean,other_c20_std = other.c20,other.c20_demean,other.c20_std

        if 'rate' in self.title and 'rate' in other.title:
            sub_c20 = self_c20 - other_c20
            sub_c20_demean = self_c20_demean - other_c20_demean
            sub_c20_std = np.sqrt(self_c20_std**2 + other_c20_std**2)

        elif ('rate' in self.title and 'rate' not in other.title) or ('rate' not in self.title and 'rate' in other.title):
            raise Exception('subtraction can not be completed between rate object and series object') 

        else:
            self_existing_solution_flag = ~self.missing_solution_flag
            other_existing_solution_flag = ~other.missing_solution_flag

            if (self_existing_solution_flag != other_existing_solution_flag).any():
                existing_solution_flag = self_existing_solution_flag & other_existing_solution_flag

                n = len(existing_solution_flag)

                assumed_self_c20 = np.zeros(n)
                assumed_other_c20 = assumed_self_c20.copy()

                assumed_self_c20[self_existing_solution_flag] = self.c20
                assumed_other_c20[other_existing_solution_flag] = other.c20

                self_c20 = assumed_self_c20[existing_solution_flag] 
                other_c20 = assumed_other_c20[existing_solution_flag] 

                assumed_self_c20_demean = np.zeros(n)
                assumed_other_c20_demean = assumed_self_c20_demean.copy()

                assumed_self_c20_demean[self_existing_solution_flag] = self.c20_demean
                assumed_other_c20_demean[other_existing_solution_flag] = other.c20_demean

                self_c20_demean = assumed_self_c20_demean[existing_solution_flag] 
                other_c20_demean = assumed_other_c20_demean[existing_solution_flag] 

                assumed_self_c20_std = np.zeros(n)
                assumed_other_c20_std = assumed_self_d20_std.copy()

                assumed_self_c20_std[self_existing_solution_flag] = self.c20_std
                assumed_other_c20_std[other_existing_solution_flag] = other.c20_std

                self_c20_std = assumed_self_c20_std[existing_solution_flag] 
                other_c20_std = assumed_other_c20_std[existing_solution_flag] 

                solution_month = self.total_month[existing_solution_flag]
                solution_counts = len(solution_month)
                missing_solution_flag = ~existing_solution_flag
                missing_month = self.total_month[missing_solution_flag]
                missing_month_counts = len(missing_month)

                info['solution_month'] = solution_month
                info['solution_counts'] = solution_counts
                info['missing_month'] = missing_month
                info['missing_month_counts'] = missing_month_counts
                info['missing_solution_flag'] = missing_solution_flag

            sub_c20 = self_c20 - other_c20
            sub_c20_demean = self_c20_demean - other_c20_demean
            sub_c20_std = np.sqrt(self_c20_std**2 + other_c20_std**2)

        return SLR_C20(info,sub_c20,sub_c20_demean,sub_c20_std)
    
    def __mul__(self,other):
        pass
    
    def __truediv__(self,other):
        pass          