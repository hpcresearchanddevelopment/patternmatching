#ifndef HAVOQGT_UPDATE_EDGE_STATE_HPP_INCLUDED
#define HAVOQGT_UPDATE_EDGE_STATE_HPP_INCLUDED

#include <unordered_set>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { namespace mpi {

template<typename Visitor>
class ues_queue {

public:
  ues_queue() {}

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

// update edge state visitor class
template<typename Graph, typename VertexData>
class ues_visitor {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;

  ues_visitor() : 
    msg_type(0) {}

  ues_visitor(vertex_locator _vertex, uint8_t _msg_type = 0) :
    vertex(_vertex), 
    msg_type(_msg_type) {}

  ues_visitor(vertex_locator _vertex, vertex_locator _parent, uint8_t _msg_type) :
    vertex(_vertex),
    parent(_parent),
    msg_type(_msg_type) {}

  ~ues_visitor() {}

  template<typename AlgData> 
  bool pre_visit(AlgData& alg_data) const {
    return false; //true;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, 
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {

    int mpi_rank = havoqgt_env()->world_comm().rank();
    int mpi_size = havoqgt_env()->world_comm().size();

    auto vertex_data = std::get<1>(alg_data)[vertex];
    auto& pattern_graph = std::get<2>(alg_data);
    //auto& vertex_active = std::get<3>(alg_data); // vertex_active
    //std::get<4>(alg_data); // edge_metadata
    //std::get<5>(alg_data); // vertex_inactive_set

    if (msg_type == 0) {
      for (size_t vertex_pattern_index = 0; 
        vertex_pattern_index < pattern_graph.vertex_data.size(); 
        vertex_pattern_index++) {
        if (pattern_graph.vertex_data[vertex_pattern_index] == vertex_data) {
          break;
          return true;
        }
      } 

      // match not found
      std::get<3>(alg_data)[vertex] = false; // vertex_active 
      // add vertex to vertex_inactive_set on this rank
      auto find_vertex = std::get<5>(alg_data).find(g.locator_to_label(vertex)); // vertex_inactive_set 
      if (find_vertex == std::get<5>(alg_data).end()) {
        auto insert_status = std::get<5>(alg_data).insert(g.locator_to_label(vertex));
        if(!insert_status.second) {
          std::cerr << "Error: failed to add an element to the set." << std::endl;
          return false;
        }     
      }    

      // send "not alive" message to neighbours
      for(eitr_type eitr = g.edges_begin(vertex); 
        eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();
        ues_visitor new_visitor(neighbor, vertex, 1);
        vis_queue->queue_visitor(new_visitor); 
      }
 
      //std::get<3>(alg_data)[vertex] = 2; TODO: reconsider this, might not be needed
  
      return true;           

    } else if (msg_type == 1) {
      //update_edge_state(alg_data); // Important : this is way too expensive, binary search will not help either 
      
      // add parent to vertex_inactive_set on this rank
      auto find_vertex = std::get<5>(alg_data).find(g.locator_to_label(parent)); // vertex_inactive_set
      if (find_vertex == std::get<5>(alg_data).end()) {
        auto insert_status = std::get<5>(alg_data).insert(g.locator_to_label(parent));
        if(!insert_status.second) {
          std::cerr << "Error: failed to add an element to the set." << std::endl;
          return false;
        }
      }

      return false; //true; ?
    } else {
      return false; //true;
    }   
    return false; //true;
  }

  friend inline bool operator>(const ues_visitor& v1, const ues_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const ues_visitor& v1, const ues_visitor& v2) {
    return false;
  }

  template<typename AlgData>  
  void update_edge_state(AlgData& alg_data) const {
    int mpi_rank = havoqgt_env()->world_comm().rank();
    auto g = std::get<0>(alg_data); // graph
    for(eitr_type eitr = g->edges_begin(vertex); 
      eitr != g->edges_end(vertex); ++eitr) {
      vertex_locator neighbour = eitr.target();
      if (g->locator_to_label(neighbour) == g->locator_to_label(parent)) {
        std::get<4>(alg_data)[eitr] = 0; // edge_metadata // 0 - inactive
        break;   
      } 
    }   
  } 
 
  vertex_locator vertex; 
  vertex_locator parent;
  uint8_t msg_type; // 0 - init, 1 - not alive
};

template <typename TGraph, typename Vertex, typename VertexData,  
  typename VertexMetaData, typename PatternGraph, typename VertexActive,
  typename EdgeMetaData>
void update_edge_state(TGraph* g, VertexMetaData& vertex_metadata, 
  PatternGraph& pattern_graph, VertexActive& vertex_active, EdgeMetaData& edge_metadata) {

  typedef typename TGraph::vertex_iterator vertex_iterator;
  typedef typename TGraph::vertex_locator vertex_locator;
  typedef typename TGraph::edge_iterator edge_iterator;

  int mpi_rank = havoqgt_env()->world_comm().rank();

  std::unordered_set<Vertex> vertex_inactive_set;
  vertex_inactive_set.clear();   

  typedef ues_visitor<TGraph, VertexData> visitor_type;
  auto alg_data = std::forward_as_tuple(g, vertex_metadata, pattern_graph, vertex_active, edge_metadata, vertex_inactive_set);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);
  vq.init_visitor_traversal_new();
  MPI_Barrier(MPI_COMM_WORLD);
  
  uint64_t global_inactive_vertex_count = havoqgt::mpi::mpi_all_reduce
    (vertex_inactive_set.size(), std::plus<size_t>(), MPI_COMM_WORLD); 
  
  std::cout << "Update Edge State | MPI Rank [" << mpi_rank << 
    "] | vertex_inactive_set Size : " << vertex_inactive_set.size() << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Update Edge State | Global Inactive Vertex Set Size : "
      << global_inactive_vertex_count << std::endl;
  }

  // set edges inactive
  for (vertex_iterator vitr = g->vertices_begin(); 
    vitr != g->vertices_end(); ++vitr) {  
    vertex_locator vertex = *vitr;
    if (vertex_active[vertex]) {
      for(edge_iterator eitr = g->edges_begin(vertex); 
        eitr != g->edges_end(vertex); ++eitr) {          
        vertex_locator neighbour = eitr.target();
        auto find_vertex = vertex_inactive_set.find(g->locator_to_label(neighbour)); 
        if (find_vertex != vertex_inactive_set.end()) {
          edge_metadata[eitr] = 0;
        }   
      }  
    } 
  }
  
  for(vertex_iterator vitr = g->delegate_vertices_begin();
    vitr != g->delegate_vertices_end(); ++vitr) {
    vertex_locator vertex = *vitr;
    if (vertex_active[vertex]) {
      for(edge_iterator eitr = g->edges_begin(vertex);
        eitr != g->edges_end(vertex); ++eitr) {
        vertex_locator neighbour = eitr.target();
        auto find_vertex = vertex_inactive_set.find(g->locator_to_label(neighbour));        
        if (find_vertex != vertex_inactive_set.end()) {
          edge_metadata[eitr] = 0;
        } 
      }   
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
} 

}} //end namespace havoqgt::mpi

#endif //HAVOQGT_UPDATE_EDGE_STATE_HPP_INCLUDED
