#ifndef HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_ITERATIVE_HPP_INCLUDED
#define HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_ITERATIVE_HPP_INCLUDED

#include <unordered_set>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

#define output_result

namespace havoqgt { namespace mpi {

template<typename IntegralType>
class vertex_state {
public:
  vertex_state() :
  //is_active(false),
  vertex_pattern_index(0) {}

  //bool is_active;
  size_t vertex_pattern_index; // TODO: change type
  std::unordered_map<size_t, IntegralType> pattern_vertex_itr_count_map; // TODO: not itr_count anymore, more like true / false 
};

template<typename Visitor>
class lppm_queue {

public:
  lppm_queue() {}

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

// label propagation pattern matching visitor class
template<typename Graph, typename VertexData>
class lppm_visitor {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;

  lppm_visitor() : 
    msg_type(0) {}

  lppm_visitor(vertex_locator _vertex, uint8_t _msg_type = 0) : 
    vertex(_vertex), 
    msg_type(_msg_type) {}

  lppm_visitor(vertex_locator _vertex, vertex_locator _parent, 
    size_t _parent_pattern_index, uint8_t _msg_type) :
    vertex(_vertex),
    parent(_parent),
    parent_pattern_index(_parent_pattern_index),
    msg_type(_msg_type) {}

  ~lppm_visitor() {}

  template<typename AlgData> 
  bool pre_visit(AlgData& alg_data) const {
    if (!std::get<4>(alg_data)[vertex]) {
      return false;
    }

    int mpi_rank = havoqgt_env()->world_comm().rank();     

    auto vertex_data = std::get<0>(alg_data)[vertex];
    //auto& pattern = std::get<1>(alg_data);
    //auto& pattern_indices = std::get<2>(alg_data);
    // std::get<4>(alg_data) - vertex_active
    auto& pattern_graph = std::get<7>(alg_data);
    // std::get<8>(alg_data) - superstep
    // std::get<9>(alg_data) - initstep 
    auto g = std::get<10>(alg_data);
    // std::get<11>(alg_data) - edge_active
    // std::get<12>(alg_data) - vertex_active_edge_set
    // std::get<12>(alg_data) - edge_metadata 

    bool match_found = false;
    bool valid_parent_found = false;

    size_t vertex_pattern_index = 0; // ID of the corresponding pattern vertex

      if (vertex.is_delegate() && g->master(vertex) != mpi_rank && msg_type == 1) { 
        // a delegate but not the controller
        // the vertex_state is only maintained on the controller

        // match vertex metadata
        // Important : no need to match edge metadata of a parent 

        // TODO: avoid figuring out vertex_pattern_index everytime
        // does vertex_data match any entry in the query pattern
        for (vertex_pattern_index = 0;
          vertex_pattern_index < pattern_graph.vertex_data.size();
          vertex_pattern_index++) {
          if (pattern_graph.vertex_data[vertex_pattern_index] == vertex_data) {
            match_found = true;   

            // verify if heard from a valid parent
            if (msg_type == 1 && match_found) {
              //match_found = false;
              for (auto e = pattern_graph.vertices[vertex_pattern_index];
                e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {
                if (pattern_graph.edges[e] == parent_pattern_index) {
                  //match_found = true;
                  valid_parent_found = true; 
                  break; 
                }
              } // for
 
              if (!valid_parent_found) {
                return false; 
              }    

            } // if            
 
          } // if

          if (valid_parent_found) {
            break;
          }
        } // for 

        // I think it never gets here 
        // initial case - return true to handle delegates // TODO: try this one too
        //if (std::get<4>(alg_data)[vertex] && msg_type == 0 && !match_found) {
        //  return true;  
        //}

        if (!match_found) {
          std::get<4>(alg_data)[vertex] = false; 
          //return true; // send to the controller ?
          return false;
        } else {
          // add parent to vertex_active_edge_set
          auto find_edge = std::get<12>(alg_data)[vertex].find(g->locator_to_label(parent));
          if (find_edge == std::get<12>(alg_data)[vertex].end()) {
            auto insert_status = std::get<12>(alg_data)[vertex].insert(g->locator_to_label(parent));
            if(!insert_status.second) {
              std::cerr << "Error: failed to add an element to the set." << std::endl;
              return false;
            }
          } else {
            // TODO: debug, why it gets here sometimes 
            //std::cerr << "(Delegate) vertex " << g->locator_to_label(vertex) 
            //  << " parent " << g->locator_to_label(parent) << " parent_pattern_index " 
            //  << parent_pattern_index << std::endl; // Test
            //std::cerr << "Error: unexpected item in the set." << std::endl;
//--            return false;
            return true;
          }
          return true; // send to the controller
        } 
 
      } 

      // for local vertex and controller only
 
      // first LP superstep of the first iteration 
      if (std::get<8>(alg_data) == 0  && std::get<9>(alg_data)) {
        for (vertex_pattern_index = 0; 
          vertex_pattern_index < pattern_graph.vertex_data.size(); 
          vertex_pattern_index++) { 
          if (pattern_graph.vertex_data[vertex_pattern_index] == vertex_data) {
             match_found = true;
             break; // Important 
          }
        }

        if (!match_found) { // TODO: controller return true?
          std::get<4>(alg_data)[vertex] = false;
          if (vertex.is_delegate() && g->master(vertex) == mpi_rank) { // controller
            //return true; // to invalidate the delegates
            return false;
          }           
          return false; 
        }
      }
     
      // if the vertex is not in the global map after the first LP 
      // superstep of the first iteration, ignore it
      if (std::get<8>(alg_data) > 0 || !std::get<9>(alg_data)) { 
        auto find_vertex = std::get<6>(alg_data).find(g->locator_to_label(vertex));
        if (find_vertex == std::get<6>(alg_data).end()) {
          return false;
        } else {
          vertex_pattern_index = find_vertex->second.vertex_pattern_index;
        }  
      }
 
      if (msg_type == 1) {
        verify_and_update_vertex_state(alg_data, vertex_pattern_index);
      }
   
      return false; 
  }
  
  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, 
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {    
    if (!std::get<4>(alg_data)[vertex]) {
      return false;
    }

    int mpi_rank = havoqgt_env()->world_comm().rank(); 

    // Important : skip this verification for the delegates as the 
    // vertex_state is only maintained on the contrller
    if (!(vertex.is_delegate() && g.master(vertex) != mpi_rank)) {
      // if the vertex is not in the global map after the first LP superstep 
      // of the first iteration, ignore it
      if (std::get<8>(alg_data) > 0 || !std::get<9>(alg_data)) {
        auto find_vertex = std::get<6>(alg_data).find(g.locator_to_label(vertex));
        if (find_vertex == std::get<6>(alg_data).end()) {
          return false; 
        }
      } 
    }

    auto vertex_data = std::get<0>(alg_data)[vertex];
    //auto& pattern = std::get<1>(alg_data);
    //auto& pattern_indices = std::get<2>(alg_data);
    // std::get<4>(alg_data) - vertex_active
    auto& pattern_graph = std::get<7>(alg_data);
    // std::get<8>(alg_data) - superstep
    // std::get<9>(alg_data) - initstep
    // std::get<10>(alg_data) - g
    // std::get<11>(alg_data) - edge_active
    // std::get<12>(alg_data) - vertex_active_edge_set    
    // std::get<13>(alg_data) - edge_metadata  

    // does vertex_data match an entry in the query pattern
    bool match_found = false;

    // TODO: do you want to compute this every time or store in the memory? Overhead is not noticable though.
    //std::vector<size_t> vertex_pattern_indices(0); // a vertex label could be a match for multiple pattern labels

    // match vertex metadata

    size_t vertex_pattern_index = 0;
    for (vertex_pattern_index = 0; 
      vertex_pattern_index < pattern_graph.vertex_data.size(); 
      vertex_pattern_index++) { 
      if (pattern_graph.vertex_data[vertex_pattern_index] == vertex_data) {
//        vertex_pattern_indices.push_back(vertex_pattern_index);
        // TODO: compare with the entry in pattern_indices to detect loop or 
        // use token passing
        match_found = true;
        break; 
      }       
    }

    if (!match_found) {
      std::get<4>(alg_data)[vertex] = false;
      return false;
      //return true; // TODO: ask Roger?
    } 

    if (msg_type == 0 && match_found) {
      // send to all the neighbours
      for(eitr_type eitr = g.edges_begin(vertex); 
        eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();      
 
        // Important : first LP superstep of the first iteration - send on all edges
        if (!(std::get<8>(alg_data) == 0  && std::get<9>(alg_data))) {
 
          if(!std::get<11>(alg_data)[eitr]) { // is edge active ?
            continue;
          }
 
        }

        // match edge metadata
        // TODO: Right now, a vertex cannot have repeating metadata. 
        // Two different vertices can have same edge metadata.
        // In order to support repeating edge metadata for the same vertex
        // edge ID need to be included in the message.  
        match_found = false;
        for (auto e = pattern_graph.vertices[vertex_pattern_index]; 
          e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {
          if (std::get<13>(alg_data)[eitr] == pattern_graph.edge_data[e]) { // edge_metadata
            match_found = true;
            break;
          }
        }

        if (!match_found) {
          std::get<11>(alg_data)[eitr] = false;
          continue;     
        }
  
        // TODO: only handling undirected grpahs, directed graphs?

//        for (auto vertex_pattern_index : vertex_pattern_indices) {
          // do this for all the pattern indices for this vertex

          lppm_visitor new_visitor(neighbor, vertex, vertex_pattern_index, 1);
          vis_queue->queue_visitor(new_visitor);

//        } // for
        
      } // for
      return true;
    } else if (msg_type == 1 && match_found) {        
      // must go all the way to the controller 
      //return true; // false?
      return false; // TODO: ask Roger?
    } else {
      return false; 
    } 
    return true;
  }

  friend inline bool operator>(const lppm_visitor& v1, const lppm_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const lppm_visitor& v1, const lppm_visitor& v2) {
    return false;
  }

  template<typename AlgData>
  uint64_t verify_and_update_vertex_state(AlgData& alg_data, 
    size_t vertex_pattern_index) const {

    typedef vertex_state<uint8_t> VertexState; // TODO: use Vertex type

    //auto& pattern = std::get<1>(alg_data);  
    //auto& pattern_indices = std::get<2>(alg_data);
    auto& pattern_graph = std::get<7>(alg_data); 
    auto g = std::get<10>(alg_data);
 
    bool match_found = false;

    // verify if parent_pattern_index is valid   
    for (auto e = pattern_graph.vertices[vertex_pattern_index]; 
      e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {  
      if (pattern_graph.edges[e] == parent_pattern_index) {
        match_found = true;
        break; 
      }  
    }

    if (!match_found) {
      return 0;
    }

    // vertex heard from a valid neighbour 
    // create an entry for this vertex in the vertex_state_map or 
    // update, if exists already 
    auto find_vertex = std::get<6>(alg_data).find(g->locator_to_label(vertex));
    if (find_vertex == std::get<6>(alg_data).end()) {
      auto insert_status = std::get<6>(alg_data).insert({g->locator_to_label(vertex), VertexState()});
      if(!insert_status.second) {
        std::cerr << "Error: failed to add an element to the map." << std::endl;
        return 0;
      }     	
      find_vertex = insert_status.first;
      find_vertex->second.vertex_pattern_index = vertex_pattern_index; // ID of the vertex in the pattern_graph 
    }

    if (std::get<6>(alg_data).size() < 1) {
      return 0;
    }

    // figure out what pattern indices are expected and add them to pattern_vertex_itr_count_map
    if (find_vertex->second.pattern_vertex_itr_count_map.size() < 1) {
      for (auto e = pattern_graph.vertices[vertex_pattern_index];
        e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {
              
        auto pattern_index = pattern_graph.edges[e];   
        
          auto find_pattern_vertex =  find_vertex->second.pattern_vertex_itr_count_map.find(pattern_index);
          if (find_pattern_vertex == find_vertex->second.pattern_vertex_itr_count_map.end()) {
            auto insert_status = find_vertex->second.pattern_vertex_itr_count_map.insert({pattern_index, 0});
            if(!insert_status.second) {
              std::cerr << "Error: failed to add an element to the map." << std::endl;
              return 0;
            }
          } 

      } // for
      
    } // if

    if (find_vertex->second.pattern_vertex_itr_count_map.size() < 1) {
      return 0;
    }

    // set status of parent_pattern_index to 1
    auto find_pattern_vertex = find_vertex->second.pattern_vertex_itr_count_map.find(parent_pattern_index);  
    if (find_pattern_vertex == find_vertex->second.pattern_vertex_itr_count_map.end()) {
      std::cerr << "Error: did not find the expected item in the map." << std::endl;
      return 0;
    }      
   
    // update status of the pattern vertex 
    if (find_pattern_vertex->second < 1) {
      find_pattern_vertex->second = 1;
    }  
   
    // add parent to vertex_active_edge_set
    auto find_edge = std::get<12>(alg_data)[vertex].find(g->locator_to_label(parent));
    if (find_edge == std::get<12>(alg_data)[vertex].end()) {
      auto insert_status = std::get<12>(alg_data)[vertex].insert(g->locator_to_label(parent));
      if(!insert_status.second) {
        std::cerr << "Error: failed to add an element to the set." << std::endl;
        return false;
      }
    } else {
      //std::string vertex_type_name = "(Local)"; // Test
      //if ((vertex.is_delegate())) { // Test
      //  vertex_type_name = "(Controller)";    
      //} 
      //std::cerr << vertex_type_name << "vertex " << g->locator_to_label(vertex) 
      //  << " parent " << g->locator_to_label(parent) << " parent_pattern_index "
      //  << parent_pattern_index << std::endl; // Test      
      //std::cerr << "Error: unexpected item in the set." << std::endl;
//--       return false;
      return 1;
    }

    return 1; 
  } 

  vertex_locator vertex;
  vertex_locator parent;
  size_t parent_pattern_index; // TODO: pass type as template argument
  uint8_t msg_type; // 0 - init, 1 - alive
};

template <typename TGraph, typename AlgData, typename VertexStateMap, 
  typename PatternGraph, typename VertexActive, typename VertexIteration,
  typename VertexSetCollection>
void verify_and_update_vertex_state_map(TGraph* g, AlgData& alg_data, 
  VertexStateMap& vertex_state_map, PatternGraph& pattern_graph, 
  VertexActive& vertex_active, 
  VertexIteration& vertex_iteration, uint64_t superstep, bool initstep, bool& global_not_finished, 
  VertexSetCollection& vertex_active_edge_set) {

  typedef typename TGraph::vertex_iterator vertex_iterator;
  typedef typename TGraph::vertex_locator vertex_locator;

  int mpi_rank = havoqgt_env()->world_comm().rank();

  // Important : invalidate vertices that have valid labels but were not added 
  // to the vertex_state_map
  if (superstep == 0 && initstep) { // Important
    for (vertex_iterator vitr = g->vertices_begin(); 
      vitr != g->vertices_end(); ++vitr) {  
      vertex_locator vertex = *vitr; 
      if (vertex_active[vertex]) {
        auto find_vertex = vertex_state_map.find(g->locator_to_label(vertex));
        if (find_vertex == vertex_state_map.end()) { 
          vertex_active[vertex] = false;    
          vertex_active_edge_set[vertex].clear();
        } 
      }  
    }

    for(vertex_iterator vitr = g->delegate_vertices_begin();
      vitr != g->delegate_vertices_end(); ++vitr) {
      vertex_locator vertex = *vitr;
      if (vertex.is_delegate() && (g->master(vertex) == mpi_rank)) {
        auto find_vertex = vertex_state_map.find(g->locator_to_label(vertex));
        if (find_vertex == vertex_state_map.end()) {
          vertex_active[vertex] = false;
        }
      }
      // skip the delegates, reduction on vertex_active will take care of them 
    }
  }

  //auto vertex_temp = vertex_state_map.begin()->first;
  //std::vector<decltype(vertex_temp)> vertex_remove_from_map_list(0); // hack
  std::vector<uint64_t> vertex_remove_from_map_list; // TODO: use Vertex type

  for (auto& v : vertex_state_map) { // TODO: use C++11 approach to remove item from map
    auto v_locator = g->label_to_locator(v.first);
    
    for (auto& p : v.second.pattern_vertex_itr_count_map) {
      if (p.second < 1) {
        vertex_remove_from_map_list.push_back(v.first);
        vertex_active[v_locator] = false; 
        break;
      } else {
        p.second = 0; // reset for next iteration   
      }   
    }   
  } // for

  if (vertex_remove_from_map_list.size() > 0) {
    global_not_finished = true;
  }

  for (auto v : vertex_remove_from_map_list) {
    if (vertex_state_map.erase(v) < 1) {
      std::cerr << "Error: failed to remove an element from the map." 
        << std::endl;  
    }    
  }

  vertex_active.all_min_reduce(); 
  MPI_Barrier(MPI_COMM_WORLD); 
}   

template <typename TGraph, typename AlgData, typename VertexActive,
  typename EdgeActive, typename VertexSetCollection>
void verify_and_update_edge_state(TGraph* g, AlgData& alg_data, 
  VertexActive& vertex_active, EdgeActive& edge_active, 
  VertexSetCollection& vertex_active_edge_set) {

  typedef typename TGraph::vertex_iterator vertex_iterator;
  typedef typename TGraph::vertex_locator vertex_locator;
  typedef typename TGraph::edge_iterator edge_iterator;

  int mpi_rank = havoqgt_env()->world_comm().rank();

  for (vertex_iterator vitr = g->vertices_begin();
    vitr != g->vertices_end(); ++vitr) {
    vertex_locator vertex = *vitr;
    if (vertex_active[vertex]) {
      for(edge_iterator eitr = g->edges_begin(vertex);
        eitr != g->edges_end(vertex); ++eitr) {
        vertex_locator neighbour = eitr.target();
        auto find_edge = vertex_active_edge_set[vertex].find(g->locator_to_label(neighbour));
        if (find_edge != vertex_active_edge_set[vertex].end()) {
          edge_active[eitr] = 1;
        } else {
          edge_active[eitr] = 0; // TODO: alternatively you could do it inside the visitor fucntion 
        }
      }
      vertex_active_edge_set[vertex].clear(); // reset 
    }
  }

  for(vertex_iterator vitr = g->delegate_vertices_begin();
    vitr != g->delegate_vertices_end(); ++vitr) {
    vertex_locator vertex = *vitr;
    if (vertex_active[vertex]) {
      for(edge_iterator eitr = g->edges_begin(vertex);
        eitr != g->edges_end(vertex); ++eitr) {
        vertex_locator neighbour = eitr.target();
        auto find_edge = vertex_active_edge_set[vertex].find(g->locator_to_label(neighbour));
        if (find_edge != vertex_active_edge_set[vertex].end()) {
          edge_active[eitr] = 1;
        } else {
          edge_active[eitr] = 0;
        } 
      }
      vertex_active_edge_set[vertex].clear(); // reset
    }
  }
  
  //vertex_active_edge_set.clear();
  MPI_Barrier(MPI_COMM_WORLD);
}

template <typename TGraph, typename VertexMetaData, typename VertexData, typename PatternData, 
  typename PatternIndices, typename VertexRank, typename VertexActive, 
  typename VertexIteration, typename VertexStateMap, typename PatternGraph, typename EdgeActive, 
  typename VertexSetCollection, typename EdgeMetadata>
void label_propagation_pattern_matching_bsp(TGraph* g, VertexMetaData& vertex_metadata, 
  PatternData& pattern, PatternIndices& pattern_indices, VertexRank& vertex_rank,
  VertexActive& vertex_active, VertexIteration& vertex_iteration, VertexStateMap& vertex_state_map, 
  PatternGraph& pattern_graph, bool initstep, bool& global_not_finished, size_t global_itr_count, 
  std::ofstream& superstep_result_file, std::ofstream& active_vertices_count_result_file, 
  EdgeActive& edge_active, EdgeMetadata& edge_metadata) {

  typedef uint64_t Vertex;

  int mpi_rank = havoqgt_env()->world_comm().rank();
  uint64_t superstep_var = 0;
  uint64_t& superstep_ref = superstep_var; 

  VertexSetCollection vertex_active_edge_set(*g);

  typedef lppm_visitor<TGraph, VertexData> visitor_type;
  auto alg_data = std::forward_as_tuple(vertex_metadata, pattern, pattern_indices, vertex_rank,
    vertex_active, vertex_iteration, vertex_state_map, pattern_graph, superstep_var, initstep, g, edge_active, 
    vertex_active_edge_set, edge_metadata);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);

  // beiginning of BSP execution
  // TODO: change for loop to use a local termination detection at the end of a seperstep
  //bool not_finished = false;
 
  for (uint64_t superstep = 0; superstep < pattern_graph.diameter; superstep++) {
    superstep_ref = superstep;
    if (mpi_rank == 0) { 
      //std::cout << "Superstep #" << superstep << std::endl;
      std::cout << "Label Propagation | Superstep #" << superstep;
    }

    double time_start = MPI_Wtime();
    //vertex_active_edge_set.clear();   
    //MPI_Barrier(MPI_COMM_WORLD); 
    
    //double time_start = MPI_Wtime();
    vq.init_visitor_traversal_new(); 
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {    
      //std::cout << "Superstep #" << superstep <<  " Synchronizing ... " << std::endl;
      std::cout <<  " | Synchronizing ...";
    }

    //vertex_active.all_min_reduce(); // do not need this here
    ///MPI_Barrier(MPI_COMM_WORLD);
 
    verify_and_update_vertex_state_map(g, alg_data, vertex_state_map, pattern_graph, 
      vertex_active, vertex_iteration, superstep, initstep, global_not_finished, 
      vertex_active_edge_set);
    //MPI_Barrier(MPI_COMM_WORLD);

    verify_and_update_edge_state(g, alg_data, vertex_active, edge_active, 
      vertex_active_edge_set);
    
    double time_end = MPI_Wtime();
    if (mpi_rank == 0) {
      //std::cout << "Superstep #" << superstep <<  " Time " << time_end - time_start << std::endl;
      std::cout << " | Time : " << time_end - time_start << std::endl;
    }

#ifdef output_result
    // result
    if (mpi_rank == 0) { 
      superstep_result_file << global_itr_count << ", LP, "
        << superstep << ", "
        << time_end - time_start << "\n"; 
    }

    // Important : This may slow things down -only for presenting results
    uint64_t active_vertices_count = 0; 
    for (auto& v : vertex_state_map) {
      auto v_locator = g->label_to_locator(v.first);
      if (v_locator.is_delegate() && (g->master(v_locator) == mpi_rank)) {
        active_vertices_count++;  
      } else if (!v_locator.is_delegate()) {
        active_vertices_count++;
      }
    }
 
    active_vertices_count_result_file << global_itr_count << ", LP, "
      << superstep << ", " 
      << active_vertices_count << "\n";
#endif

    // TODO: global reduction on global_not_finished before next iteration

  } // for 
  // end of BSP execution  
}  

}} //end namespace havoqgt::mpi

#endif //HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_ITERATIVE_HPP_INCLUDED 
