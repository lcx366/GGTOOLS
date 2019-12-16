class Nodes(object):
    '''
    class Nodes
    - attributes:
        - nodes
        - nodes_index
        - region
    
    - methods:
    '''
    def __init__(self,nodes,nodes_index,region):
        self.nodes = nodes
        self.nodes_index = nodes_index
        self.region = region
        
    def __repr__(self):
        return 'Nodes for mascons'