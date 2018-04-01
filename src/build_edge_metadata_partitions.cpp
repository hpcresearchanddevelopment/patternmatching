#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/interprocess/managed_heap_memory.hpp>

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/environment.hpp>

#include <havoqgt/edge_data_db.hpp>

namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;
using namespace havoqgt;

typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;

template<typename T>
  using SegmentAllocator = bip::allocator<T, segment_manager_t>;

template <typename DelegateGraph, typename VertexLocator, typename EdgeIterator>
void populate_edge_metadata(DelegateGraph* graph, VertexLocator& vertex, 
  std::string vertex_class, std::ifstream& original_edge_metadata_file,
  std::ofstream& edge_metadata_file) {
  for(EdgeIterator eitr = graph->edges_begin(vertex);
    eitr != graph->edges_end(vertex); ++eitr) {
    VertexLocator neighbour = eitr.target();
    uint64_t edge_metadata = 555;
    edge_metadata_file << vertex_class 
    << "," << graph->locator_to_label(vertex)
    << "," << graph->locator_to_label(neighbour)
    << "," << edge_metadata << "\n";
  }
}

int main(int argc, char** argv) {
  typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;

  int mpi_rank(0), mpi_size(0);

  // havoqgt_init
  havoqgt::havoqgt_init(&argc, &argv);
  {

  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  havoqgt::get_environment();

  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    havoqgt::get_environment().print();
    //print_system_info(false);
  }
  MPI_Barrier(MPI_COMM_WORLD);  

  std::string graph_input = argv[1];
  std::string backup_filename = argv[2];
  std::string metadata_input_dir = argv[3];
  std::string metadata_output_dir = argv[4];
  std::string vertex_metadata_input_filename = argv[5];
  std::string result_dir = argv[6];
   
  // TODO: parse commandline

  MPI_Barrier(MPI_COMM_WORLD);  

  // graph
  if(backup_filename.size() > 0) {
    distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
  }

  havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());

  //segment_manager_t* segment_manager = ddb.get_segment_manager();
  //  bip::allocator<void, segment_manager_t> alloc_inst(segment_manager);

  if (mpi_rank == 0) {
    std::cout << "Loading Graph ... " << std::endl;
  }

  //graph_type *graph = segment_manager->
  //  find<graph_type>("graph_obj").first;
  //assert(graph != nullptr);

  graph_type *graph = ddb.get_segment_manager()->
    find<graph_type>("graph_obj").first;
  assert(graph != nullptr);

  // edge data
  if (mpi_rank == 0) {
    std::cout << "Loading / Initializing Edge Data ... " << std::endl;
  }

  typedef uint8_t edge_data_type; 

  // TODO: figure out a way to get it from graph_type
  // see edge_data_value_type in parallel_edge_list_reader.hpp

  typedef graph_type::edge_data<edge_data_type, 
    bip::allocator<edge_data_type, segment_manager_t>> edge_data_t;

  edge_data_t* edge_data_ptr = ddb.get_segment_manager()->
    find<edge_data_t>("graph_edge_data_obj").first;
  //assert(edge_data_ptr != nullptr);  
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Done Loading Graph." << std::endl;
  }

  //graph->print_graph_statistics(); // causes MPI error
  MPI_Barrier(MPI_COMM_WORLD);

  // application parameters
  bool do_output_vertex_data = false;  

  /////////////////////////////////////////////////////////////////////////////
  
  // build edge metadata partitions 
  {
    //if(mpi_rank == 0) {
    //  std::cout << "Building Edge Metadata Partitions." << std::endl;
    //}
 
    typedef typename graph_type::vertex_locator vertex_locator; 
    typedef typename graph_type::vertex_iterator vertex_iterator;
    typedef typename graph_type::edge_iterator edge_iterator;

    typedef uint64_t Vertex;
    typedef uint64_t Edge;
    typedef uint64_t VertexData;  
    typedef uint64_t EdgeData;

    typedef graph_type::vertex_data<VertexData, std::allocator<VertexData> > VertexMetaData;
    typedef graph_type::edge_data<EdgeData, std::allocator<EdgeData> > EdgeMetadata;

    double time_start = MPI_Wtime();
    double time_end = MPI_Wtime();  

    // vertex containers
    VertexMetaData vertex_metadata(*graph);
    
    // edge containers
    EdgeMetadata edge_metadata(*graph); 

    if(mpi_rank == 0) { 
      std::cout << "Application | Allocated Miscellaneous Containers" << std::endl;
    }
    
    // application parameters // TODO: commandline input
    std::string edge_metadata_input_filename = metadata_input_dir; // TODO: add support fo multiple input files

    MPI_Barrier(MPI_COMM_WORLD);
 
    /////////////////////////////////////////////////////////////////////////////
    
    // build the distributed vertex data db
    time_start = MPI_Wtime();

    //if (use_degree_as_vertex_data) {
    //  vertex_data_db_degree<graph_type, VertexMetaData, Vertex, VertexData>
    //    (graph, vertex_metadata);   
    //} else { 
      edge_data_db<graph_type, Vertex, VertexData, EdgeData, VertexMetaData, EdgeMetadata>
        (graph, vertex_metadata, vertex_metadata_input_filename, 
        edge_metadata, edge_metadata_input_filename, 10000);
//        edge_data_db<graph_type, Vertex, VertexData, VertexMetaData>
//          (graph, vertex_metadata, vertex_metadata_input_filename, 10000); 
        // TODO: each rank reads 10K lines from file at a time
    //} 

    MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this?
    time_end = MPI_Wtime();
    if(mpi_rank == 0) {
      std::cout << "Application Time | Vertex Data DB : " 
        << time_end - time_start << std::endl;
    }

    if (do_output_vertex_data) {
      std::string vertex_data_filename = result_dir + "/" +
      std::to_string(0) + "/all_ranks_vertex_data/vertex_data_" + std::to_string(mpi_rank);
      std::ofstream vertex_data_file(vertex_data_filename, std::ofstream::out);

      for (vertex_iterator vitr = graph->vertices_begin(); vitr != graph->vertices_end();
        ++vitr) {
        vertex_locator vertex = *vitr;
        vertex_data_file << mpi_rank << ", l, " << graph->locator_to_label(vertex) 
          << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";
        if(mpi_rank == 0) { // Test  
          std::cout << vertex_metadata[vertex] << std::endl;
        } 
      } 	
    	
      for(vertex_iterator vitr = graph->delegate_vertices_begin();
        vitr != graph->delegate_vertices_end(); ++vitr) {
        vertex_locator vertex = *vitr;

        if (vertex.is_delegate() && (graph->master(vertex) == mpi_rank)) {
          vertex_data_file << mpi_rank << ", c, " << graph->locator_to_label(vertex) 
            << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";
          if(mpi_rank == 0) { // Test
            std::cout << vertex_metadata[vertex] << std::endl;
          }
        } else {	
          vertex_data_file << mpi_rank << ", d, " << graph->locator_to_label(vertex) 
            << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";  
          if(mpi_rank == 0) { // Test
            std::cout << vertex_metadata[vertex] << std::endl;
          }
        }
      }
  
      vertex_data_file.close();
    } 

    /////////////////////////////////////////////////////////////////////////////

    bool do_write_edgelist_partitions = false ;

    if (do_write_edgelist_partitions) {

    std::string edgelist_filename = metadata_output_dir + 
      "/all_ranks_edgelist/edgelist_" + std::to_string(mpi_rank);
    std::ofstream edgelist_file(edgelist_filename, std::ofstream::out);

    // local vertices
    for (vertex_iterator vitr = graph->vertices_begin();
      vitr != graph->vertices_end(); ++vitr) {          
      vertex_locator vertex = *vitr;      
      //edgelist_file << mpi_rank << ", l, " << graph->locator_to_label(vertex) 
      //  << ", " << graph->degree(vertex) << "\n";    
      for(edge_iterator eitr = graph->edges_begin(vertex);
        eitr != graph->edges_end(vertex); ++eitr) {
        vertex_locator neighbour = eitr.target();
        edgelist_file << " l," << graph->locator_to_label(vertex)
        << "," << graph->locator_to_label(neighbour) << "\n";    
      }      
    }

    // delegate vertices
    for(vertex_iterator vitr = graph->delegate_vertices_begin();
      vitr != graph->delegate_vertices_end(); ++vitr) {
      vertex_locator vertex = *vitr;
      std::string vertex_type_name = "d"; 	      
      if (vertex.is_delegate() && (graph->master(vertex) == mpi_rank)) {
        //edgelist_file << mpi_rank << ", c, " << graph->locator_to_label(vertex)
        //  << ", " << graph->degree(vertex) << "\n";
        vertex_type_name = "c";
      } //else { 
        //edgelist_file << mpi_rank << ", d, " << graph->locator_to_label(vertex)
        //  << ", " << graph->degree(vertex) << "\n";    
      //} 
      for(edge_iterator eitr = graph->edges_begin(vertex);
        eitr != graph->edges_end(vertex); ++eitr) {
        vertex_locator neighbour = eitr.target();
        edgelist_file << vertex_type_name << "," << graph->locator_to_label(vertex)
        << "," << graph->locator_to_label(neighbour) << "\n";
      } 
    }  

    edgelist_file.close();

    } // do_write_edgelist_partitions

    else {

    if(mpi_rank == 0) {
      std::cout << "Building Edge Metadata Partitions." << std::endl;
    }

    std::string original_edge_metadata_filename = metadata_input_dir; // TODO: add support fo multiple input files
    std::ifstream original_edge_metadata_file(original_edge_metadata_filename, 
      std::ifstream::in); 

    std::string edgelist_filename = metadata_output_dir +
      "/all_ranks_edgelist/edgelist_" + std::to_string(mpi_rank);
    std::ifstream edgelist_file(edgelist_filename, std::ifstream::in);
 
    std::string edge_metadata_filename = metadata_output_dir +
    "/all_ranks_edge_metadata/edge_metadata_" + std::to_string(mpi_rank);
    std::ofstream edge_metadata_file(edge_metadata_filename, std::ofstream::out);

    // metadata - read original, write for partitions
//    std::string line;
//    while (std::getline(original_edge_metadata_file, line)) {
          
//    }

    std::string vertex_class = "l";
    // local vertices
    for (vertex_iterator vitr = graph->vertices_begin();
      vitr != graph->vertices_end(); ++vitr) {          
      vertex_locator vertex = *vitr;      
      //edge_metadata_file << mpi_rank << ", l, " << graph->locator_to_label(vertex) 
      //  << ", " << graph->degree(vertex) << "\n";    
//      for(edge_iterator eitr = graph->edges_begin(vertex);
//        eitr != graph->edges_end(vertex); ++eitr) {
//        vertex_locator neighbour = eitr.target();
//        uint64_t edge_metadata = 555;
//        edge_metadata_file << " l," << graph->locator_to_label(vertex)
//        << "," << graph->locator_to_label(neighbour) 
//        << "," << edge_metadata << "\n";    
//      }      

      populate_edge_metadata<graph_type, vertex_locator, edge_iterator>
        (graph, vertex, vertex_class, original_edge_metadata_file, 
        edge_metadata_file); 
    }

    // delegate vertices
    for(vertex_iterator vitr = graph->delegate_vertices_begin();
      vitr != graph->delegate_vertices_end(); ++vitr) {
      vertex_locator vertex = *vitr;
//      std::string vertex_type_name = "d"; 	      
      vertex_class = "d";
      if (vertex.is_delegate() && (graph->master(vertex) == mpi_rank)) {
        //edge_metadata_file << mpi_rank << ", c, " << graph->locator_to_label(vertex)
        //  << ", " << graph->degree(vertex) << "\n";
//        vertex_type_name = "c";
        vertex_class = "c";
      } //else { 
        //edge_metadata_file << mpi_rank << ", d, " << graph->locator_to_label(vertex)
        //  << ", " << graph->degree(vertex) << "\n";    
      //} 
//      for(edge_iterator eitr = graph->edges_begin(vertex);
//        eitr != graph->edges_end(vertex); ++eitr) {
//        vertex_locator neighbour = eitr.target();
//        uint64_t edge_metadata = 555; 
//        edge_metadata_file << vertex_type_name << "," << graph->locator_to_label(vertex)
//        << "," << graph->locator_to_label(neighbour) 
//        << "," << edge_metadata << "\n";
//      } 
      populate_edge_metadata<graph_type, vertex_locator, edge_iterator>
        (graph, vertex, vertex_class, original_edge_metadata_file,
        edge_metadata_file);
    } 

    edge_metadata_file.close();
    edgelist_file.close();
    original_edge_metadata_file.close();
 
    } // else
    
  } // build edge metadata partitions

  if(mpi_rank == 0) {
    std::cout << "Done." << std::endl;
  }

  } // havoqgt_init
  // END Main MPI
  havoqgt::havoqgt_finalize();

  return 0;  
} // main
