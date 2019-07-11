

class disjoint_set:
    "nodes are mapped by their id to a tuple of their key and rank"

    def __init__(self, ids = None):
        "initializes the disjoint set. ids is an optional list of ids to add to the set"
        self.node_map = {}
        if ids != None:
            for id in ids:
                self.add(id)

    def add(self, id):
        "adds an id to the disjoint set"
        if not id in self.node_map:
            self.node_map[id] = (id,0)
    
    def find(self, id):
        "finds the key of a given value in the disjoint set"
        if not id in self.node_map:
            raise ValueError("the id {id} does not exist in the disjoint_set")

        #see if node points to itself
        key = self.node_map[id][0]
        if key == id:
            return self.node_map[key]
        else:
            new_tuple = self.find(key)
            new_key = new_tuple[0]
            new_rank = new_tuple[1]

            #apply path compression
            self.node_map[key] = (new_key, new_rank)
            return self.node_map[key]

    def getKey(self, id):
        "gets the key of the corresponding id"

        if id not in self.node_map:
            raise ValueError("the id {id} does not exist")
        return self.find(id)[0]

    def union(self, left_id, right_id):
        "merges two disjoint sets"

        if(not left_id in self.node_map):
            raise ValueError("the left id {left_id} does not exist")
        if(not right_id in self.node_map):
            raise ValueError("the left id {right_id} does not exist")

        #the two ids are already in the same set
        if(self.getKey(left_id) == self.getKey(right_id)):
            return

        current_left_rank = self.node_map[left_id][1]
        left_tuple = self.find(left_id)
        left_key = left_tuple[0]
        left_rank = left_tuple[1]

        current_right_rank = self.node_map[right_id][1]
        right_tuple = self.find(right_id)
        right_key = right_tuple[0]
        right_rank = right_tuple[1]

        if left_rank == right_rank:
            self.node_map[right_id] = (left_key, current_right_rank)
            self.node_map[left_id] = (left_key, current_left_rank)

            self.node_map[right_key] = (left_key, right_rank)
            self.node_map[left_key] = ((left_key, left_rank + 1))
        elif left_rank>right_rank:
            self.node_map[right_id] = (left_key, current_right_rank)
            self.node_map[left_id] = (left_key, current_left_rank)

            self.node_map[right_key] = (left_key, right_rank)
            self.node_map[left_key] = ((left_key, left_rank + 1))
        else:
            self.node_map[left_id] = (right_key, current_left_rank)
            self.node_map[right_id] = (right_key, current_right_rank)

            self.node_map[right_key] = (right_key, right_rank)
            self.node_map[left_key] = ((right_key, left_rank + 1))

    def equivalent(self, left_id, right_id):
        "determines if two ids are in the same disjoint set"

        if(not left_id in self.node_map):
            return False
        if(not right_id in self.node_map):
            return False
        return self.getKey(left_id) == self.getKey(right_id)

