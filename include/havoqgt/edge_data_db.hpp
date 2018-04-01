#ifndef HAVOQGT_EDGE_DATA_DB_HPP_INCLUDED
#define HAVOQGT_EDGE_DATA_DB_HPP_INCLUDED

//#define BOOST_FILESYSTEM_VERSION 3
//#define BOOST_FILESYSTEM_NO_DEPRECATED 

#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { namespace mpi {

template<typename Visitor>
class vertex_data_queue {

public:
  vertex_data_queue() {}

  bool push(Visitor const& element) {
    data.push_back(element);
    return true;
  }

  void pop() {
    data.pop_back();
  }
 
  Visitor const& top() {
    return data.back();
  } 
  
  size_t size() const {
    return data.size();;
  }

  bool empty() const {
    return data.empty();
  }

  void clear() {
    data.clear();
  }

protected:
  std::vector<Visitor> data;
};

template<typename Graph, typename VertexData, typename EdgeData>
class vertex_data_visitor {

public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;

  vertex_data_visitor() : 
    vertex_data(0),
    do_update_vertex_data(false) {}

  vertex_data_visitor(vertex_locator _vertex) : 
    vertex(_vertex), 
    vertex_data(0), 
    do_update_vertex_data(false) {}

  vertex_data_visitor(vertex_locator _vertex, VertexData _vertex_data) :
    vertex(_vertex), 
    vertex_data(_vertex_data), 
    do_update_vertex_data(true) {}

  vertex_data_visitor(vertex_locator _vertex, vertex_locator _parent, EdgeData _edge_data) :
    vertex(_vertex),
    parent(_parent),
    edge_data(_edge_data),
    do_update_vertex_data(true) {}
 
  ~vertex_data_visitor() {}

  template<typename AlgData> 
  bool pre_visit(AlgData& alg_data) const {
    return true;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, 
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  /*template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    if (!do_update_vertex_data) {
      //std::cout << "False" << std::endl; // test 
      return false;
    } 

    std::get<0>(alg_data)[vertex] = vertex_data;    
    //std::cout << "Visiting " << g.locator_to_label(vertex) << " " << std::get<0>(alg_data)[vertex] << std::endl; // test  
    return true;
  }*/
 
  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    if (!do_update_vertex_data) {
      //std::cout << "False" << std::endl; // test 
      return false;
    }

    int mpi_rank = havoqgt_env()->world_comm().rank();
    bool match_found = false;

    //if (vertex.is_delegate()) {
    //  return false; 
    //}      

    //std::get<0>(alg_data)[vertex] = edge_data;    
    //std::cout << "Visiting " << g.locator_to_label(vertex) << " " << std::get<0>(alg_data)[vertex] << std::endl; // test

    if ((g.locator_to_label(vertex) == 5002504) && mpi_rank == 10) {
   
    for(eitr_type eitr = g.edges_begin(vertex); 
      eitr != g.edges_end(vertex); ++eitr) {
      vertex_locator neighbor = eitr.target();
      //std::cout << mpi_rank << " " << g.locator_to_label(vertex) << " " << g.locator_to_label(neighbor) << std::endl; // test
      if (g.locator_to_label(neighbor) == 4808458 ) { //g.locator_to_label(parent)) {
        std::get<0>(alg_data)[eitr] = edge_data;  
        match_found = true;
        // update edge metadata here
        //break;
      }      
    }

    // binary search 4808458

    }

    //std::binary_search(g.edges_begin(vertex), g.edges_end(vertex), g.local_degree(vertex));     
    
//    if (g.local_degree(vertex) < 1) {
//      return false;
//    } else if (g.local_degree(vertex) == 1) {
      // update edge metadata  
//      return true;
//    }


    //uint64_t low = 0;
    //uint64_t high = g.local_degree(vertex) - 1;
    //uint64_t mid = (low + high) / 2;

    //eitr_type eitr_high = g.edges_begin(vertex) + high; // Test
    //std::cout << "N " << g.locator_to_label(eitr_high.target()) << std::endl; // Test
    //std::cout << "P " << g.locator_to_label(parent) << std::endl; // Test

    //eitr_type eitr_low = g.edges_begin(vertex);
    //eitr_type eitr_high = g.edges_begin(vertex) + (g.local_degree(vertex) - 1);

    //if (eitr_high == g.edges_end(vertex)) { // Test // works
    //  std::cout << "T" << std::endl;
    //} else {
    //  std::cout << "False" << std::endl;
    //}
      
/*    while(low <= high) {
      eitr_type eitr = g.edges_begin(vertex) + mid;  

      //if (eitr == g.edges_end(vertex)) {
      //  std::cout << "END" << std::endl;  
      //  break; 
      //}
 
      vertex_locator neighbor = eitr.target();

      if (g.locator_to_label(neighbor) == g.locator_to_label(parent)) {
        match_found = true;
        // update edge metadata here
        break;    
      } else if (g.locator_to_label(neighbor) < g.locator_to_label(parent)) {
        low = mid + 1;
      } else {
        if (mid == 0) {
          //std::cout << "MID = 0" << std::endl;
          break;  
        }    
        high = mid - 1;           
      } 
      mid = (low + high) / 2; 
    } // while*/

    if (match_found) {
      std::cout << " >>> Found" << std::endl;
    }

    if (vertex.is_delegate() && !match_found) {
      return true;    
    } else if (!vertex.is_delegate() && !match_found) {
      //std::cout << low << " " << high << " " << mid << " " << "Not found!" << std::endl;  
      //std::cerr << "Error: target for the local vertex not found." << std::endl;
      //return false;
    } else {
      return true;
    }
    
    return true;
  }
 
  friend inline bool operator>(const vertex_data_visitor& v1, const vertex_data_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const vertex_data_visitor& v1, const vertex_data_visitor& v2) {
    return false;
  }

  //bool vertex_comparator(Graph& g, vertex_locator& a, vertex_locator& b) {
  //  return g.locator_to_label(a) < g.locator_to_label(b);
  //}

  vertex_locator vertex;
  vertex_locator parent;
  VertexData vertex_data;
  EdgeData edge_data;
  bool do_update_vertex_data;
};

template <typename TGraph, typename Vertex, typename VertexData,
  typename EdgeData, typename VertexMetaData, typename EdgeMetadata, 
  typename VertexEntry>
void edge_data_db(TGraph* g, VertexMetaData& vertex_metadata, 
  EdgeMetadata& edge_metadata, VertexEntry& vertex_entry) {
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  typedef vertex_data_visitor<TGraph, VertexData, EdgeData> visitor_type;
  auto alg_data = std::forward_as_tuple(edge_metadata);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);

  bool visitor_created = false; // false;

//  std::cout << "MPI Rank " << mpi_rank << " queuing local vertex entries ... " << std::endl;  
  for (auto entry : vertex_entry) {
    auto source = std::get<0>(entry);
    auto target = std::get<1>(entry);
    auto edge_data = std::get<2>(entry);
    auto source_location = g->label_to_locator(source);   
    auto target_location = g->label_to_locator(target);	

    //if (mpi_rank == 0) { // Test
    //  std::cout
    //  << source << " " << target << " " << edge_data << std::endl;
    //}

    //visitor_type new_visitor(source_location, edge_data);
    if (source == 5002504 && !visitor_created && mpi_rank == 10) {
    visitor_type source_visitor(source_location, target_location, edge_data);
    vq.queue_visitor(source_visitor); 
    visitor_created = true; 
    }
//    visitor_type target_visitor(target_location, source_location, edge_data);
//    vq.queue_visitor(target_visitor); 
  }
//  std::cout << "MPI Rank " << mpi_rank << " done queuing ... building distributed vertex data db ... " << std::endl;
  
  // TODO: implement a more efficient 'visitor_traversal'; e.g., only visit the ones already in the queue
  //MPI_Barrier(MPI_COMM_WORLD); // TODO: deadlock
  vq.init_visitor_traversal_new();
  MPI_Barrier(MPI_COMM_WORLD);
//  std::cout << "MPI Rank " << mpi_rank << " is done." << std::endl; 
}

bool get_files_in_dir(boost::filesystem::path& dir_path, 
  std::vector<std::string>& file_paths, std::string wildcard) {
  const boost::regex filename_filter(wildcard + ".*");  

  if(!boost::filesystem::exists(dir_path) || 
    !boost::filesystem::is_directory(dir_path)) {
    std::cerr << "Error: Invalid directory path." << std::endl;
    return false;    
  } else {
    boost::smatch what; 
    boost::filesystem::directory_iterator end_itr;
    for (boost::filesystem::directory_iterator itr(dir_path);
      itr != end_itr; ++itr) {
      if (!boost::filesystem::is_regular_file(itr->status())) {
        continue; 
      } else if ( !boost::regex_match(itr->path().filename().string(), what, filename_filter) ) {
        continue;
      } else { 
        file_paths.push_back(itr->path().string()); 
      }
    }
     
    if (file_paths.size() < 1) { 
      return false; 
    } else {
      return true;
    } 
  }  
} 

template <typename TGraph, typename Vertex, typename VertexData,
  typename EdgeData, typename VertexMetaData, typename EdgeMetadata>
void read_file_and_build_edge_data_db(std::string vertex_data_input_filename, 
  size_t chunk_size, TGraph* g, VertexMetaData& vertex_metadata, 
  EdgeMetadata& edge_metadata) {

  int mpi_rank = havoqgt_env()->world_comm().rank();

  typedef std::vector<std::tuple<Vertex, Vertex, VertexData>> VertexEntry;
  VertexEntry vertex_entry;

  std::ifstream vertex_data_input_file(vertex_data_input_filename,
  std::ifstream::in);
  std::string line;
  while (std::getline(vertex_data_input_file, line)) {
    std::istringstream iss(line);
    //std::string source_vertex_class_name = "";
    Vertex source(0), target(0);    
    VertexData edge_data(0);
    //iss >> source >> target >> edge_data; // replaced with boost functions

    std::vector<std::string> tokens(0);  
    boost::split(tokens, line, boost::is_any_of(", "), boost::token_compress_on);
    
    source = boost::lexical_cast<uint64_t>(tokens[1].c_str()); 
    target = boost::lexical_cast<uint64_t>(tokens[2].c_str());
    edge_data = boost::lexical_cast<uint64_t>(tokens[3].c_str());    

    //if (mpi_rank == 0) { // Test
    //  std::cout 
    //  << source << " " << target << " " << edge_data << std::endl;
    //} 
    
    vertex_entry.push_back(std::forward_as_tuple(source, target, edge_data));
  }
  vertex_data_input_file.close();

  edge_data_db<TGraph, Vertex, VertexData, EdgeData, VertexMetaData, EdgeMetadata, VertexEntry>
    (g, vertex_metadata, edge_metadata, vertex_entry);

  // TODO: process in chunks. 
  // Important: same issue as number of files per-rank.
  // For now, make sure each vertex data file is not too big.
  // It can handle any number of files, so reading files in chunks is not 
  // important. 
}

template <typename TGraph, typename Vertex, typename VertexData, 
  typename EdgeData, typename VertexMetaData, typename EdgeMetadata>
void edge_data_db(TGraph* g, VertexMetaData& vertex_metadata, 
  std::string base_filename, EdgeMetadata& edge_metadata, 
  std::string edge_metadata_base_filename, size_t chunk_size) {
//template <typename TGraph, typename Vertex, typename VertexData,typename VertexMetaData>
//void edge_data_db(TGraph* g, VertexMetaData& vertex_metadata,
//  std::string base_filename, size_t chunk_size) {

  int mpi_rank = havoqgt_env()->world_comm().rank();
  int mpi_size = havoqgt_env()->world_comm().size();

  if (mpi_rank == 0) {
    std::cout << "Building distributed vertex data db ... " << std::endl;
  }

  boost::filesystem::path dir_path(base_filename);
  std::string wildcard = dir_path.filename().string();
  dir_path.remove_filename();
 
  std::vector<std::string> file_paths;  

  if(!get_files_in_dir(dir_path, file_paths, wildcard)) {    
    std::cerr << "Error: Failed to read input files." << std::endl;    
    return; 
  }

  // Important: if not all the ranks have the same number of files to process 
  size_t max_files_per_rank = 0;
  if (file_paths.size() < mpi_size) {
    max_files_per_rank = 1;  
  } else {
      if (file_paths.size() % mpi_size == 0) {
      max_files_per_rank = file_paths.size() / mpi_size; 
    } else {
      double quotient = file_paths.size() / (double) mpi_size; 
      max_files_per_rank = static_cast<size_t>(floor(quotient) + 1);
    }
  }
  if (mpi_rank == 0) {
    std::cout << "Maximum number of files per-rank " << 
      max_files_per_rank << std::endl; 
  }

  for (size_t i = mpi_rank <= file_paths.size() - 1 ? mpi_rank : 
    static_cast<size_t>(mpi_rank % file_paths.size()), j = 0; 
    //i < file_paths.size(); i+=mpi_size) {
    j < max_files_per_rank; j++) {

    assert(i >= 0 && i < file_paths.size());    
//    std::cout << "MPI Rank " << mpi_rank << " processing file [" << i << "] " 
//      << file_paths[i] << " ... " << std::endl;

    read_file_and_build_edge_data_db<TGraph, Vertex, VertexData, EdgeData, VertexMetaData, EdgeMetadata>
      (file_paths[i], chunk_size, g, vertex_metadata, edge_metadata);

    if (i + mpi_size < file_paths.size()) {
      i+=mpi_size;  
    }     
  } 

  if (mpi_rank == 0) {
    std::cout << "Done building vertex data db." << std::endl;
  }   
    
} 

}} //end namespace havoqgt::mpi
 

#endif //HAVOQGT_EDGE_DATA_DB_HPP_INCLUDED
