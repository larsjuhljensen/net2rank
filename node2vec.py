'''
    This script is used to run the node2vec algorithm on the STRING functional networks.
'''

import argparse
import os
import warnings
import numba
from gensim.models import Word2Vec
from gensim.models.callbacks import CallbackAny2Vec
from pecanpy.cli import pecanpy
import numpy as np
import h5py
from typing import List, Tuple, Iterable, Union
import gzip
from multiprocessing import Pool
import itertools

class H5pyData:

    @staticmethod
    def write(proteins: Union[np.ndarray, List, Tuple],
            embedding: np.ndarray,
            save_path: str,
            precision: int,
            chunk_size: int = 10000) -> None:
        '''
        Write proteins and embeddings to HDF5 file efficiently.
        
        Args:
            proteins: Array-like of protein identifiers
            embedding: numpy array of embeddings (n_proteins x embedding_dim)
            save_path: Path to save the HDF5 file
            precision: Precision of embeddings (16 or 32)
            chunk_size: Size of chunks for HDF5 dataset
        '''
        # Convert proteins to numpy array if needed
        proteins = np.array(proteins).astype('U').reshape(-1)
        embedding = np.array(embedding)
        
        # Validate inputs
        if len(proteins) != len(embedding):
            raise ValueError(f"Number of proteins ({len(proteins)}) doesn't match number of embeddings ({len(embedding)})")
        
        if precision not in [16, 32]:
            raise ValueError(f"Precision must be 16 or 32, got {precision}")
        
        # Determine dtype for embeddings
        dtype = np.float16 if precision == 16 else np.float32
        embedding = embedding.astype(dtype)
        
        n_proteins = len(proteins)
        embedding_dim = embedding.shape[1]
        
        with h5py.File(save_path, 'w') as f:
            # Create groups to organize data
            f.create_group('metadata')
            
            # Store metadata
            f['metadata'].attrs['n_proteins'] = n_proteins
            f['metadata'].attrs['embedding_dim'] = embedding_dim
            f['metadata'].attrs['precision'] = precision
            
            # Create datasets with chunking and compression
            protein_ds = f.create_dataset(
                'proteins',
                shape=(n_proteins,),
                dtype=h5py.string_dtype(),
                chunks=(min(chunk_size, n_proteins),),
                compression='gzip',
                compression_opts=4
            )
            
            embedding_ds = f.create_dataset(
                'embeddings',
                shape=(n_proteins, embedding_dim),
                dtype=dtype,
                chunks=(min(chunk_size, n_proteins), embedding_dim),
                compression='gzip',
                compression_opts=4
            )
            
            # Write data in chunks for better memory management
            for i in range(0, n_proteins, chunk_size):
                end_idx = min(i + chunk_size, n_proteins)
                
                protein_chunk = proteins[i:end_idx]
                embedding_chunk = embedding[i:end_idx]
                
                protein_ds[i:end_idx] = protein_chunk
                embedding_ds[i:end_idx] = embedding_chunk

    @staticmethod
    def read(file_path: str,precision: int) -> tuple[np.ndarray, np.ndarray]:
        '''
        Read proteins and embeddings from HDF5 file.
        
        Args:
            file_path: Path to HDF5 file
            precision: Precision of embeddings (16 or 32)
            
        Returns:
            tuple: (proteins array, embeddings array)
        '''
        with h5py.File(file_path, 'r') as f:
            proteins = f['proteins'][:]
            proteins = np.vectorize(lambda x: str(x)[2:-1])(proteins)
            embeddings = f['embeddings'][:]
        if precision == 16:
            embeddings = embeddings.astype(np.float16)
        elif precision == 32:
            embeddings = embeddings.astype(np.float32)

        return proteins, embeddings
    

    
    

def query_single_species(querys, precision, aligned_path):

    # logger.info(f'Loading {aligned_path}...')
    proteins,embeddings = H5pyData.read(aligned_path,precision)

    protein2index = {protein: index for index, protein in enumerate(proteins)}

    indices = []
    for query in querys:
        try:
            index = protein2index[query]
            indices.append( index)
        except NameError as e:
            logger.info(f'Query {query} not found. Error: {e}')
    
    ## get the embeddings
    query_proteins = [proteins[index] for index in indices]
    query_embeddings = [embeddings[index] for index in indices]

    return query_proteins, query_embeddings

def query_embedding(querys:Iterable, query_dir:str,precision:int, n_jobs:int) -> Tuple[np.ndarray,np.ndarray]:
    '''
        Query the embeddings for the given proteins in multiple species.
        Args:
            querys: Iterable, the proteins to query.
            aligned_dir: str, the directory of the aligned embeddings.
            precision: int, the precision of the embeddings, either 32 or 16.
            n_jobs: int, the number of jobs to run in parallel.
        Returns:
            output_proteins: np.array, the proteins that are found.
            output_embeddings: np.array, the embeddings of the proteins.

    '''

    ## check if the directory of query_dir exists
    if not os.path.exists(query_dir):
        raise FileNotFoundError(f'{query_dir} does not exist')
    
    if precision not in [32,16]:
        raise ValueError('Precision should be either 32 or 16')

    ## put the quries in a dict where the key is taxid and the value is the list of querys
    query_dict = {}

    for query in querys:
        taxid = query.split('.')[0]
        if taxid not in query_dict:
            query_dict[taxid] = [query]
        else:
            query_dict[taxid].append(query)

    ## for each taxid get the data and the querys
    with Pool(n_jobs) as p:
        results = p.starmap(query_single_species, [(query_dict[taxid],precision, f'{query_dir}/{taxid}.h5') for taxid in query_dict])

    output_proteins = list()
    output_embeddings = list()

    for result in results:
        output_proteins.append(result[0])
        output_embeddings.append(result[1])

    output_proteins = np.array(list(itertools.chain.from_iterable(output_proteins)))
    output_embeddings = np.array(list(itertools.chain.from_iterable(output_embeddings)))

    return output_proteins, output_embeddings


## a class to handle the gz file

class GzipData:
    @staticmethod
    def string2idx(file_path:str,temp_path)->dict:

        nodes = dict()

        ## check if the directory of temp_path exists
        if not os.path.exists(os.path.dirname(temp_path)):
            os.makedirs(os.path.dirname(temp_path))

        edges_writer = csv.writer(open(temp_path, 'w'), delimiter='\t')

        with gzip.open(file_path, 'rt') as f:
            
            reader = csv.reader(f, delimiter=' ')

            next(reader) # skip the header

            for row in reader:
                
                if row[0] not in nodes:
                    nodes[row[0]] = len(nodes)
                if row[1] not in nodes:
                    nodes[row[1]] = len(nodes)
                
                ## get the index of the protein
                src_idx = nodes[row[0]]
                dst_idx = nodes[row[1]]
                weight = int(row[-1])/1000

                edges_writer.writerow([src_idx,dst_idx,weight])

        return nodes
    
    @staticmethod
    def read_nodes(file_path:str)->dict:
        nodes = dict()
        with gzip.open(file_path, 'rt') as f:
            reader = csv.reader(f, delimiter=' ')
            
            next(reader) # skip the header

            for row in reader:
                if row[0] not in nodes:
                    nodes[row[0]] = len(nodes)
                if row[1] not in nodes:
                    nodes[row[1]] = len(nodes)
        return nodes


class callback(CallbackAny2Vec):
    '''Callback to print loss after each epoch.'''

    def __init__(self,total_epoch):
        self.epoch = 1
        self.loss_to_be_subed = 0
        self.total_epoch = total_epoch
        self.saved_loss = list()

    def on_epoch_end(self, model):
        loss = model.get_latest_training_loss()
        loss_now = loss - self.loss_to_be_subed
        self.loss_to_be_subed = loss
        self.saved_loss.append([self.epoch,loss_now])
        self.epoch += 1



def learn_embeddings(walks,epochs=1,dimensions=64,window_size=10,workers=1,random_state=1234,taxid="",writer=None,hs=0):
    """_summary_

    Args:
        epochs (int, optional): _description_. Defaults to 1.
        walks (int, optional): _description_. Defaults to 20.
        dimensions (int, optional): _description_. Defaults to 64.
        window_size (int, optional): _description_. Defaults to 10.
        workers (int, optional): _description_. Defaults to 1.
        random_state (int, optional): _description_. Defaults to 1234.
        output_path (str, optional): _description_. Defaults to ''.

    Returns:
        _type_: _description_
    """    
    # cb = callback(epochs,taxid,writer) ## de-comment this line if you want to save the loss

    model = Word2Vec(
        walks,
        vector_size=dimensions,
        window=window_size,
        min_count=0,
        sg=1,
        hs=hs,
        workers=workers,
        epochs=epochs,
        seed=random_state,
        compute_loss=True,
        # callbacks=[cb]  ## de-comment this line if you want to save the loss
    )

    return model





def check_mode(g, mode,weighted,p,q):
    """Check mode selection.

    Give recommendation to user for pecanpy mode based on graph size and density.

    """
    # mode = args.mode
    # weighted = args.weighted
    # p = args.p
    # q = args.q

    # Check unweighted first order random walk usage
    if mode == "FirstOrderUnweighted":
        if not p == q == 1 or weighted:
            raise ValueError(
                f"FirstOrderUnweighted only works when weighted = False and "
                f"p = q = 1, got {weighted=}, {p=}, {q=}",
            )
        return

    if mode != "FirstOrderUnweighted" and p == q == 1 and not weighted:
        warnings.warn(
            "When p = 1 and q = 1 with unweighted graph, it is highly "
            f"recommended to use FirstOrderUnweighted over {mode} (current "
            "selection). The runtime could be improved greatly with improved  "
            "memory usage.",
        )
        return

    # Check first order random walk usage
    if mode == "PreCompFirstOrder":
        if not p == q == 1:
            raise ValueError(
                f"PreCompFirstOrder only works when p = q = 1, got {p=}, {q=}",
            )
        return

    if mode != "PreCompFirstOrder" and p == 1 == q:
        warnings.warn(
            "When p = 1 and q = 1, it is highly recommended to use "
            f"PreCompFirstOrder over {mode} (current selection). The runtime "
            "could be improved greatly with low memory usage.",
        )
        return

    # Check network density and recommend appropriate mode
    g_size = g.num_nodes
    g_dens = g.density
    if (g_dens >= 0.2) & (mode != "DenseOTF"):
        warnings.warn(
            f"Network density = {g_dens:.3f} (> 0.2), it is recommended to use "
            f"DenseOTF over {mode} (current selection)",
        )
    if (g_dens < 0.001) & (g_size < 10000) & (mode != "PreComp"):
        warnings.warn(
            f"Network density = {g_dens:.2e} (< 0.001) with {g_size} nodes "
            f"(< 10000), it is recommended to use PreComp over {mode} (current "
            "selection)",
        )
    if (g_dens >= 0.001) & (g_dens < 0.2) & (mode != "SparseOTF"):
        warnings.warn(
            f"Network density = {g_dens:.3f}, it is recommended to use "
            f"SparseOTF over {mode} (current selection)",
        )
    if (g_dens < 0.001) & (g_size >= 10000) & (mode != "SparseOTF"):
        warnings.warn(
            f"Network density = {g_dens:.3f} (< 0.001) with {g_size} nodes "
            f"(>= 10000), it is recommended to use SparseOTF over {mode} "
            "(current selection)",
        )



def read_graph(path,p,q,workers,verbose,weighted,directed,extend,gamma,random_state,mode,delimiter,implicit_ids):
    """Read input network to memory.

    Depending on the mode selected, reads the network either in CSR
    representation (``PreComp`` and ``SparseOTF``) or 2d numpy array
    (``DenseOTF``).

    """

    if directed and extend:
        raise NotImplementedError("Node2vec+ not implemented for directed graph yet.")

    if extend and not weighted:
        print("NOTE: node2vec+ is equivalent to node2vec for unweighted graphs.")

    # if task in ["tocsr", "todense"]:  # perform conversion then save and exit
    #     g = graph.SparseGraph() if task == "tocsr" else graph.DenseGraph()
    #     g.read_edg(path, weighted, directed, delimiter)
    #     g.save(output)
    #     exit()

    pecanpy_mode = getattr(pecanpy, mode, None)
    g = pecanpy_mode(p, q, workers, verbose, extend, gamma, random_state)

    if path.endswith(".npz"):
        g.read_npz(path, weighted, implicit_ids=implicit_ids)
    else:
        g.read_edg(path, weighted, directed, delimiter)

    check_mode(g, mode,weighted,p,q)

    return g


def preprocess(g):
    """Preprocessing transition probabilities with timer."""
    g.preprocess_transition_probs()


def simulate_walks(num_walks, walk_length, g):
    """Simulate random walks with timer."""
    return g.simulate_walks(num_walks, walk_length)


class PecanpyEmbedder():

    def __init__(self,graph_path,p=1,q=1,workers=-1,
                 weighted=True,directed=False,
                 extend=False,gamma=0,random_state=1234,
                 delimiter:str='\t'):
        super().__init__()

        if workers == -1:
            workers = numba.config.NUMBA_DEFAULT_NUM_THREADS

        ## load the graph
        self.graph = read_graph(graph_path,p,q,workers,verbose=False,
                                weighted=weighted,directed=directed,
                                extend=extend,gamma=gamma,
                                random_state=random_state,
                                mode="SparseOTF",delimiter=delimiter,
                                implicit_ids=False)
        preprocess(self.graph)


    def generate_walks(self,num_walks:int,walk_length:int) -> list:
        return simulate_walks(num_walks,walk_length,self.graph)
        



    def learn_embeddings(self, walks, epochs=1,dimensions=128, 
                         window_size=5,workers=-1,
                         negative=5,
                         hs=0,sg=1,
                         random_state=1234) -> Word2Vec:
        """
            Word2Vec API of gensim

            Parameters
            ----------
            walks : list, list of walks
            epochs : int, number of epochs, default 1.
            dimensions : int, number of dimensions, default 128.
            window_size : int, window size, default 5.
            workers : int, number of workers, default -1 (all workers).
            negative : int, if >0, negative sampling will be used, number of negative samples if use negative sampling, default 5.
            hs : int, if 1, use hierarchical softmax; if 0, use negative sampling, default 0.
            sg : int, if 1, use skip-gram; if 0, use CBOW, default 1.
            random_state : int, random state, default 1234.
        """
       
        cb = callback(epochs)

        if workers == -1:
            workers = numba.config.NUMBA_DEFAULT_NUM_THREADS

        model = Word2Vec(
            walks,
            vector_size=dimensions,
            window=window_size,
            min_count=0,
            sg=sg,
            hs=hs,
            workers=workers,
            epochs=epochs,
            seed=random_state,
            compute_loss=True,
            callbacks=[cb],
            negative=negative
        )

        return model



def run_single_embedding(species_file, temp_path, output_folder, dimensions, 
                         p, q, num_walks, walk_length, window_size, sg, 
                         negative, epochs, workers, random_state):
    """Run single species-specific embedding."""
    ## process the gz file
    logger.info(f"Processing {species_file}...")
    nodes = GzipData.string2idx(species_file,temp_path)


    if len(nodes) > 50000:
        logger.warning(f"Number of nodes in {species_file} is {len(nodes)}, if it fails, try larger memory")

    logger.info(f"Running embedding for {species_file}...")
    # Read the graph
    embedder = PecanpyEmbedder(temp_path,p=p,q=q,workers=workers,weighted=True,directed=False,
                               extend=False,gamma=0,random_state=random_state,delimiter='\t')
    # Generate the walks

    embeddings = embedder.learn_embeddings(embedder.generate_walks(num_walks=num_walks,walk_length=walk_length),
                                      epochs=epochs,dimensions=dimensions,window_size=window_size,
                                      workers=workers,negative=negative,hs=0,sg=sg,random_state=random_state)
    
    emb = embeddings.wv.vectors
    index = embeddings.wv.index_to_key

    proteins = list(nodes.keys())

    ## map the index to the protein
    map_proteins = [proteins[int(i)] for i in index]

    ## save the embeddings
    species = species_file.split('/')[-1].split('.')[0]
    save_path = f"{output_folder}/{species}.h5"
    H5pyData.write(map_proteins,emb,save_path,32)
    logger.info(f"Embedding for {species_file} is saved at {save_path}")

    return None

def main(args):

    species_name = os.path.basename(args.input_network).split(".")[0]

    # create a temporary file
    temp_path = f"{args.node2vec_output}/{species_name}.tsv"

    # run the single embedding
    run_single_embedding(args.input_network, temp_path, args.node2vec_output, args.dimensions, 
                        args.p, args.q, args.num_walks, args.walk_length, args.window_size, args.sg, 
                        args.negative, args.epochs, args.workers, args.random_state)

    # remove the temporary file
    os.remove(temp_path)

    return None


if __name__ == "__main__":


    parser_single_embedding = argparse.ArgumentParser(description="Run node2vec on STRING functional networks with PecanPy.")

    parser_single_embedding.add_argument('-i','--input_network', type=str, help='File to run the\
                                        embedding for, e.g. <species_name>.tsv.gz. During running, the file will be processed into another temporary tsv file,\n\
                                        proteins will be replaced by integers, and the scores will be converted to float (0~1). Once finished, the temporary file will be deleted.')
                                      
    parser_single_embedding.add_argument('-o','--node2vec_output', type=str, help='Path to the output folder to save the embeddings.\n\
                                         The embeddings will be saved in the format: <output_folder>/<species_name>.h5')
    
    
    ### model parameters, optional
    parser_single_embedding.add_argument('-d', '--dimensions', type=int, default=128,help='The number of dimensions for the embedding.')
    parser_single_embedding.add_argument('-p', '--p', type=float, default=0.3, help='The return parameter for the random walk.')
    parser_single_embedding.add_argument('-q', '--q', type=float, default=0.7, help='The in-out parameter for the random walk.')
    parser_single_embedding.add_argument('--num_walks', type=int, default=10, help='The number of walks to perform.')
    parser_single_embedding.add_argument('--walk_length', type=int, default=50,help='The length of the walk.')
    parser_single_embedding.add_argument('--window_size', type=int, default=5, help='The window size for the skip-gram model.')
    parser_single_embedding.add_argument('--sg', type=int, default=1, help='The type of training to use for the skip-gram model. 0 for cbow, 1 for skip-gram.')
    parser_single_embedding.add_argument('--negative', type=int, default=5, help='The number of negative samples to use for training the model.')
    parser_single_embedding.add_argument('-e', '--epochs', default=5, type=int, help='The number of epochs to train the model.')
    parser_single_embedding.add_argument('--workers', type=int, default=-1, help='The number of workers to use for training the model.')
    parser_single_embedding.add_argument('--random_state', type=int, default=1234, help='The random state to use for the random number generator.')

    args = parser_single_embedding.parse_args()

    main(args)