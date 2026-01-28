from abc import ABC, abstractmethod
import pandas as pd


class FileParser(ABC):
    def __init__(self, file_path):
        self.file_path = file_path

    @abstractmethod
    def parse(self):
        pass


class OBOParser(FileParser):

    def parse(self):
        rows = []
        current_term = None
        obsolete = False

        with open(self.file_path) as f:
            for line in f:
                line = line.strip()

                if line == "[Term]":
                    # save previous term
                    if current_term and not obsolete:
                        rows.append(current_term)

                    # start new term
                    current_term = {
                        "go_id": "",
                        "name": "",
                        "namespace": "",
                        "parents": [],
                        "relationships": [],   # NEW
                        "definition": None,    # NEW
                        "synonyms": []         # NEW
                    }
                    obsolete = False

                elif current_term is not None:

                    if line.startswith("id:"):
                        current_term["go_id"] = line.split("id:")[1].strip()

                    elif line.startswith("name:"):
                        current_term["name"] = line.split("name:")[1].strip()

                    elif line.startswith("namespace:"):
                        current_term["namespace"] = line.split("namespace:")[1].strip()

                    elif line.startswith("is_a:"):
                        parent_id = line.split("is_a:")[1].split()[0]
                        current_term["parents"].append(parent_id)

                    elif line.startswith("relationship:"):
                        # Example: relationship: part_of GO:0008150
                        rel = line.split("relationship:")[1].strip()
                        current_term["relationships"].append(rel)

                    elif line.startswith("def:"):
                        current_term["definition"] = line.split("def:")[1].strip()

                    elif line.startswith("synonym:"):
                        current_term["synonyms"].append(
                            line.split("synonym:")[1].split('"')[1]
                        )

                    elif line.startswith("is_obsolete: true"):
                        obsolete = True

        # save last term
        if current_term and not obsolete:
            rows.append(current_term)

        return pd.DataFrame(rows) # Create a table where each row is one GO term and each column is one attribute of the term.


class GAFParser(FileParser):

    def parse(self):
        rows = []

        with open(self.file_path) as f:
            for line in f:
                line = line.strip()

                if line.startswith("!"):
                    continue

                parts = line.split("\t")

                if len(parts) < 17:
                    continue

                current_row = {}

                if True:
                    current_row["gene_id"] = parts[1]
                    current_row["gene_symbol"] = parts[2]
                    current_row["qualifier"] = parts[3]
                    current_row["go_id"] = parts[4]
                    current_row["evidence"] = parts[6]
                    current_row["aspect"] = parts[8]
                    current_row["gene_name"] = parts[9]
                    current_row["taxon"] = parts[12].replace("taxon:", "")

                rows.append(current_row)

        return pd.DataFrame(rows)
    
class Term: #Represents a single GO term.
    def __init__(self, go_id, name, namespace, definition=None, synonyms=None, is_a=None):
        self.go_id = go_id
        self.name = name
        self.namespace = namespace
        self.is_a = is_a or [] # Raw relationships from parsing (GO IDs)
        
        self.definition = definition #idk if we'll need this
        self.synonyms = synonyms or [] #idk if we'll need this

        self._parents = set()
        self._children = set()
        
        self._relationships = []
        #_parents and _children are protected because you don’t want external code modifying them directly.
        
    def add_parent(self, parent: 'Term'):
        self._parents.add(parent)
        parent._children.add(self)
        
    def add_relationship(self, relationship: 'Relationship'):# no classes
        self._relationships.append(relationship)  


            
    def __repr__(self):
        return f"GOTerm({self.go_id}, {self.name}, {self.namespace})"
#Now we have full collection of terms => we can start interconnecting them

class TermCollection: 
    def __init__(self):
        self._terms = {}

    def add_term(self, term: Term):
        self._terms[term.go_id] = term

    def build_vertical_relationship(self): # it creates relationships from the is_a thing
            for term in self._terms.values():
                for parent_id in term.is_a:
                    parent = self.terms.get(parent_id)
                    if parent != None:
                        term.add_parent(parent)
                        
    def build_horizontal_relationships(self, regulation_data):
        pass
                    
    def get_term(self, go_id):
        return self.terms.get(go_id)
    
    def get_parents(self, go_id):
        term = self.get_term(go_id)
        return term.parents if term else []

    def get_children(self, go_id):
        term = self.get_term(go_id)
        return term.children if term else []   
    
    def get_ancestors(self, go_id): # Return a set of all ancestor terms of a given GO term (recursively)
        term = self.get_term(go_id)
        if not term:
            return set()

        ancestors = set()

        def explore(term):
            for parent in term.parents:
                if parent not in ancestors:
                    ancestors.add(parent)
                    explore(parent)  
        explore(term)
        return ancestors 
    
    def get_descendants(self, go_id): # Return a set of all descendant terms of a given GO term (recursively)
            term = self.get_term(go_id)
            if not term:
                return set()
    
            descendants = set()
    
            def explore(term):
                for child in term.children:
                    if child not in descendants:
                        descendants.add(child)
                        explore(child)
    
            explore(term)
            return descendants
        
class GeneAnnotation:
    def __init__(self, gene, go_id, qualifier=None, aspect=None, evidence=None):
        self.gene = gene
        self.go_id = go_id
        self.qualifier = qualifier
        self.aspect = aspect
        self.evidence = evidence
        # convert GAF aspect to branch
        aspect_map = {'P': BiologicalProcess, 'F': MolecularFunction, 'C': CellularComponent}
        self.branch = aspect_map.get(aspect, Branch)()
        
    def link_term(self, term):
        """Assign a Term object to this annotation after resolution."""
        term = term_collection.get_term(self.go_id)
        if term != None:
            self.term = term
            
class AnnotationCollection: # holds all annotations and provides query functions
    def __init__(self):
        # just a list of all annotations
        self.annotations = []

    def add_annotation(self, annotation):
        """Add a GeneAnnotation object to the collection."""
        self.annotations.append(annotation)

    def link_terms(self, term_collection):
        """Resolve all GO IDs to Term objects."""
        for ann in self.annotations:
            ann.link_term(term_collection)

    def get_by_gene(self, gene_id):
        """Return all annotations for a specific gene."""
        return [ann for ann in self.annotations if ann.gene.gene_id == gene_id]

    def get_by_term(self, go_id):
        """Return all annotations for a specific GO term."""
        return [ann for ann in self.annotations if ann.go_id == go_id]

    def get_by_aspect(self, aspect):
        """Return all annotations with a specific aspect (F, P, C)."""
        return [ann for ann in self.annotations if ann.aspect == aspect]

    def get_by_evidence(self, evidence):
        """Return all annotations with a specific evidence code."""
        return [ann for ann in self.annotations if ann.evidence == evidence]

    def __len__(self):
        return len(self.annotations)

    def __iter__(self):
        return iter(self.annotations)

    def __repr__(self):
        return f"<AnnotationCollection: {len(self.annotations)} annotations>"
