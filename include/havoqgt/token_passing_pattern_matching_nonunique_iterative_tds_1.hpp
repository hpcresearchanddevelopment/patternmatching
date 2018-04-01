#ifndef HAVOQGT_TOKEN_PASSING_PATTERN_MATCHING_NONUNIQUE_ITERATIVE_TDS_1_HPP_INCLUDED 
#define HAVOQGT_TOKEN_PASSING_PATTERN_MATCHING_NONUNIQUE_ITERATIVE_TDS_1_HPP_INCLUDED

#include <algorithm>
#include <array>
#include <deque>
#include <limits>
#include <unordered_set>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

// TODO: improve
static bool enable_token_agreegation = true;
static bool enable_vertex_token_source_cache = false;
static bool path_checking_filter = true;
static uint64_t visitor_count;

namespace havoqgt { namespace mpi {

template<typename Visitor>
class tppm_queue_tds {

public:
  tppm_queue_tds() {}

  bool push(Visitor const& element) {
    data.push_back(element);
    return true;
  }

  void pop() {
    //data.pop_back();
    data.pop_front();
  }

  Visitor const& top() {
    //return data.back();
    return data.front();
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
  //std::vector<Visitor> data;
  std::deque<Visitor> data;
};

// token passing pattern matching visitor class
template<typename Graph, typename Vertex, typename BitSet>
class tppm_visitor_tds {

public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;

  tppm_visitor_tds() : 
    itr_count(0), 
    do_pass_token(false),
    is_init_step(true),
    ack_success(false),
    msg_type(0), 
    //source_index_pattern_indices(0), 
    parent_pattern_index(0) {
      visited_vertices[itr_count] = vertex;  
  }

  tppm_visitor_tds(vertex_locator _vertex) :  
    vertex(_vertex), 
    itr_count(0),
    do_pass_token(false), 
    is_init_step(true),
    ack_success(false),
    msg_type(0),
    //source_index_pattern_indices(0), 
    parent_pattern_index(0) {
      visited_vertices[itr_count] = vertex;  
  }
  
  template <typename VertexLocatorArrayStatic> 
  tppm_visitor_tds(vertex_locator _vertex,
    vertex_locator _parent,  
    vertex_locator _target_vertex, 
    Vertex _vertex_label,
    Vertex _target_vertex_label, 
    uint64_t _sequence_number, 
    VertexLocatorArrayStatic& _visited_vertices, 
    size_t _itr_count, 
    //size_t _max_itr_count, 
    //size_t _source_index_pattern_indices, 
    size_t _parent_pattern_index, 
    //bool _expect_target_vertex = true, 
    bool _do_pass_token = true, 
    bool _is_init_step = false, 
    bool _ack_success = false, 
    uint8_t _msg_type = 1) :
 
    vertex(_vertex),
    parent(_parent),
    target_vertex(_target_vertex),
    vertex_label(_vertex_label),
    target_vertex_label(_target_vertex_label),
    sequence_number(_sequence_number),  
    itr_count(_itr_count), 
    //max_itr_count(_max_itr_count), 
    //expect_target_vertex(_expect_target_vertex), 
    do_pass_token(_do_pass_token), 
    is_init_step(_is_init_step),
    ack_success(_ack_success),
    msg_type(_msg_type), 
    //source_index_pattern_indices(_source_index_pattern_indices), 
    parent_pattern_index(_parent_pattern_index) {
      if (!_ack_success) {
        // copy to visited_vertices from the one received from the parent
        std::copy(std::begin(_visited_vertices), 
          std::end(_visited_vertices), // TODO: std::begin(_visited_vertices) + itr_count ? 
          std::begin(visited_vertices));  
        visited_vertices[itr_count + 1] = vertex; // Important : itr_count for the next hop            
      }
  }  

  template<typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    visitor_count++;

    if(!std::get<13>(alg_data)[vertex]) {
      return false;
    }

    int mpi_rank = havoqgt_env()->world_comm().rank();

    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto& pattern = std::get<1>(alg_data);
    auto& pattern_indices = std::get<2>(alg_data);
    //auto& pattern_graph = std::get<4>(alg_data);
    auto max_itr_count = std::get<7>(alg_data);
    auto expect_target_vertex = std::get<8>(alg_data);
    auto g = std::get<11>(alg_data); // graph
    // std::get<12>(alg_data); // vertex_token_source_set   
    // std::get<13>(alg_data); // vertex_active
    // std::get<14>(alg_data); // template_vertices
    // std::get<15>(alg_data); // vertex_active_edges_map
    // std::get<16>(alg_data); // pattern_selected_vertices
    // std::get<17>(alg_data); // paths_result_file
    // std::get<18>(alg_data); // pattern_enumeration_indices
    // std::get<19>(alg_data); // visitor_set
    auto superstep = std::get<20>(alg_data); // superstep_var
    // std::get<21>(alg_data); // vertex_sequence_number
    // std::get<22>(alg_data); // pattern_aggregation_steps 
 
    auto source_index_pattern_indices = 0;  

    auto is_vertex_delegate_slave = false;
    if (vertex.is_delegate() && g->master(vertex) != mpi_rank) {
      is_vertex_delegate_slave = true;
    }

//    if (superstep > 0) { // Test
//      std::cerr << "Superstep " << superstep
//        << " Vertex " << g->locator_to_label(vertex) << " "
//        << " is delegate slave " << is_vertex_delegate_slave    
//        << " pre-visiting ... "
//        << std::endl; // Test
//    } // Test 

    if (ack_success) {
      // TODO: does controller need to return true?
      return true; // Important : delegate forwarding to the controller  
    } 

    // TODO: remove or replace   
    if(enable_vertex_token_source_cache) { 
    // verify if this vertex has already forwarded a token from the originating vertex
    if (!is_init_step && max_itr_count > itr_count) {
      auto find_token_source_forwarded = std::get<12>(alg_data)[vertex]
        .find(g->locator_to_label(target_vertex));
      if (find_token_source_forwarded != std::get<12>(alg_data)[vertex].end()) {
        return false;
      }
    }
    } //enable_vertex_token_source_cache
    
    if (!do_pass_token && is_init_step && itr_count == 0 && superstep == 0) { 
      // probably it never gets here
//      if (!(find_vertex->second.vertex_pattern_index == pattern_indices[0] 
//        && vertex_data == pattern[0])) {
      if(vertex_data == pattern[0]) { 
        return true;
      } else {
        return false; 
      } 
    } else if (!is_init_step) { 
       // forward tokens
       
       // TODO: remove, // not needed for tds  
       // target_vertex cannot relay the token
       //if (do_pass_token && (max_itr_count > itr_count) &&
       //  (g->locator_to_label(vertex) == g->locator_to_label(target_vertex))) {
       //  return false;
       //}  

       auto new_itr_count = itr_count + 1;
       auto next_pattern_index = source_index_pattern_indices + new_itr_count; // expected next pattern_index
       // TODO: notice next_pattern_index is redundent now, same as new_itr_count; just use new_itr_count       
       auto vertex_pattern_index = 0; //find_vertex->second.vertex_pattern_index;

       // TODO: read vertex_pattern_index from template_vertices
       // vertex_data == pattern[next_pattern_index] and 
       // vertex_pattern_index == pattern_indices[next_pattern_index] must hold  
       // otherwise return false
        
       // verify vertex data  
       // verify if received from a valid parent
       if (vertex_data != pattern[next_pattern_index] && 
         parent_pattern_index != pattern_indices[next_pattern_index - 1]) {
         return false;
       }

       bool match_found = false; 
       
       BitSet vertex_template_vertices(std::get<14>(alg_data)[vertex]); // template_vertices
  
       if (vertex_template_vertices.none()) {
         return false;  
       }

       // verify template vertex match  
       // TODO: replace the loop by vertex_template_vertices.test(pattern_indices[next_pattern_index])  
       for (size_t i = 0; i < vertex_template_vertices.size(); i++) {  
         if (vertex_template_vertices.test(i)) {
           if (i == pattern_indices[next_pattern_index]) {
             match_found = true;
             vertex_pattern_index = i; 
             break;  
           }    
         }   
       }
 
       if (!match_found) {
         return false;
       }    

       // verify if received from a valid parent // TODO: remove redundent checks
       if (vertex_data == pattern[next_pattern_index] &&
         vertex_pattern_index == pattern_indices[next_pattern_index]) {
         if (vertex_data == pattern[next_pattern_index] && 
           parent_pattern_index == pattern_indices[next_pattern_index - 1]) { 

//+           if (vertex.is_delegate() && g->master(vertex) != mpi_rank) { // delegate but not the controller
//+             return true;
//+           }  

           if (do_pass_token && ( (max_itr_count > itr_count) || (max_itr_count == itr_count) ) ) { 
             // TODO: no agreegation on a delegate?   
             //if (is_vertex_delegate_slave) {
             //  return true; // Important : delegate forwarding to the controller
             //} 
 
             if (enable_token_agreegation) {
               auto find_token = std::get<19>(alg_data)->find(*this);
               if (find_token == std::get<19>(alg_data)->end()) {
                 //std::cout << "MPI Rank: " << mpi_rank << " Superstep " << superstep << " "  
                 //  << g->locator_to_label(vertex) << " " << g->locator_to_label(target_vertex)
                 //  << " Not in the visitor set." << std::endl; // Test
               } else {
//                 std::cout << "MPI Rank: " << mpi_rank << " Superstep " << superstep << " "
//                   << " is delegate slave " << is_vertex_delegate_slave << " "
//                   << g->locator_to_label(vertex) << " " << g->locator_to_label(target_vertex)
//                   << " Already in the visitor set." << std::endl; // Test  
                 return false;
               }                  
             }
             // enable_token_agreegation 
             
             // TODO: remove or replace     
             // OK to forward the token, now update vertex_token_source_set   
             if (enable_vertex_token_source_cache) { 
             auto find_token_source_forwarded = std::get<12>(alg_data)[vertex]
               .find(g->locator_to_label(target_vertex));
             if (find_token_source_forwarded == std::get<12>(alg_data)[vertex].end()) {
               auto insert_status = std::get<12>(alg_data)[vertex]
                 .insert(g->locator_to_label(target_vertex));
               if (!insert_status.second) {
                 std::cerr << "Error: failed to add an element to the set." << std::endl;
                 return false;
               }
               // std::cout << g.locator_to_label(vertex) << " adding " 
               //   << g.locator_to_label(target_vertex)
	       //  << " to the vertex set" << std::endl; // Test     
	     } else {
               std::cerr << "Error: unexpected item in the set." << std::endl;
               return false;
             }
             } // enable_vertex_token_source_cache

             // TDS 
             //std::get<18>(alg_data); // pattern_enumeration_indices
             //std::cout << g->locator_to_label(vertex) << " pattern_enumeration_indices " 
             //  << new_itr_count << ", " 
             //  << std::get<18>(alg_data)[new_itr_count] << std::endl; // Test
             if (std::get<18>(alg_data)[new_itr_count] == new_itr_count) {
               // duplicate is not expected in visited_vertices  
               for (size_t i = 0; i < new_itr_count; i++) { // Important : must be < not <=
                 if (g->locator_to_label(visited_vertices[i]) == g->locator_to_label(vertex)) {
                   //break;
                   //std::cout << g->locator_to_label(vertex) << " returning false from pre_visit ... " << std::endl; // Test 
                   return false;
                 }
               }  
             } else if (std::get<18>(alg_data)[new_itr_count] < new_itr_count) {
               // vertex is expected in visited_vertices
               if (g->locator_to_label(visited_vertices[std::get<18>(alg_data)[new_itr_count]]) 
                 != g->locator_to_label(vertex)) {
                 return false;
               }
             } else {
               std::cerr << "Error: invalid value." << std::endl;
               return false;
             }
         
             // add to the visitor set 
             { 
             //auto find_token = std::get<19>(alg_data)->find(*this);
             auto insert_status = std::get<19>(alg_data)->insert(*this);  
             if (!insert_status.second) {
               std::cerr << "MPI Rank: " << mpi_rank << " Superstep " << superstep << " "
                 << g->locator_to_label(vertex) << " " << g->locator_to_label(target_vertex)
                 << " Error: failed to add an element to the set." << std::endl; // Test
               //std::cerr << "Error: failed to add an element to the set." << std::endl;
               //return false;
             } else {
//               std::cout << "MPI Rank: " << mpi_rank << " Superstep " << superstep << " " 
//                 << g->locator_to_label(vertex) << " " << g->locator_to_label(target_vertex)
//                 << " is delegate slave " << is_vertex_delegate_slave << " "  
//                 << " Added to the visitor set." 
                 //<< " Pointer address " << std::get<19>(alg_data) 
//                 << " size " << std::get<19>(alg_data)->size()   
//                 << std::endl; // Test
             }             
             }                       

           } // max_itr_count > itr_count || max_itr_count == itr_count 

//           if (do_pass_token && (max_itr_count == itr_count)) {
              
           /*if (enable_token_agreegation) {
             //if (g->locator_to_label(vertex) == g->locator_to_label(target_vertex)) {
             //  std::cout << "MPI Rank: " << mpi_rank << " BC " << std::endl; // Test
             //}

             auto find_token = std::get<19>(alg_data)->find(*this);

             //if (g->locator_to_label(vertex) == g->locator_to_label(target_vertex)) {
             //  std::cout << "MPI Rank: " << mpi_rank << " AD " << std::endl; // Test
             //}

             if (find_token == std::get<19>(alg_data)->end()) {
               std::cout << "MPI Rank: " << mpi_rank << " "
                 << g->locator_to_label(vertex) << " " << g->locator_to_label(target_vertex)
                 << " Not in the visitor set." << std::endl; // Test
             } else {
               std::cout << "MPI Rank: " << mpi_rank << " " 
                 << g->locator_to_label(vertex) << " " << g->locator_to_label(target_vertex)
                 << " Already in the visitor set." << std::endl; // Test  
               //return false;
             }                  
           } // enable_token_agreegation
              
           { 
             //auto find_token = std::get<19>(alg_data)->find(*this);
             auto insert_status = std::get<19>(alg_data)->insert(*this);  
             if (!insert_status.second) {
               std::cerr << "MPI Rank: " << mpi_rank << " "
                 << g->locator_to_label(vertex) << " " << g->locator_to_label(target_vertex)
                 << " Error: failed to add an element to the set." << std::endl; // Test
               //std::cerr << "Error: failed to add an element to the set." << std::endl;
               //return false;
             } else {
               std::cout << "MPI Rank: " << mpi_rank << " "
                 << g->locator_to_label(vertex) << " " << g->locator_to_label(target_vertex) 
                 << " Added to the visitor set." << std::endl; // Test
             }             
           }*/

             // TODO: 
             // local vertex and controller process the visitor here
             // no need to add to the visitor_set
             // return false
             // delegate forwards to the controller
//             return true;  
//           } // max_itr_count == itr_count                   

           if (is_vertex_delegate_slave) {                         
             return true; // Important : delegate forwarding to the controller 
           } else if (!is_vertex_delegate_slave && (max_itr_count == itr_count)) { 
             return true; // Important  
           } else {
             return false; // Important
             //return true; // Test 
           }

         } else {
           return false;
         }    
       } else {
         return false;
       }  
    } else {
      return false;  
    }
//++    return true;
    return false;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
    AlgData& alg_data) const {
    if(!std::get<13>(alg_data)[vertex]) {
      return false;
    } else { 
      return visit(g, vis_queue, alg_data);
    }
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    //visitor_count++;
    
    if(!std::get<13>(alg_data)[vertex]) {
      return false;
    }

    int mpi_rank = havoqgt_env()->world_comm().rank();
    
    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto& pattern = std::get<1>(alg_data);
    auto& pattern_indices = std::get<2>(alg_data);   
    //auto& pattern_graph = std::get<4>(alg_data);

    //auto pattern_cycle_length = std::get<7>(alg_data);   
    auto max_itr_count = std::get<7>(alg_data);
    auto expect_target_vertex = std::get<8>(alg_data);
    //auto pattern_valid_cycle = std::get<8>(alg_data);
    //auto& pattern_found = std::get<9>(alg_data);
    //auto& edge_metadata = std::get<10>(alg_data); 
    // std::get<11>(alg_data) // graph
    // std::get<12>(alg_data); // vertex_token_source_set
    // std::get<13>(alg_data); // vertex_active 
    // std::get<14>(alg_data); // template_vertices
    // std::get<15>(alg_data); // vertex_active_edges_map
    // std::get<16>(alg_data); // pattern_selected_vertices
    // std::get<17>(alg_data); // paths_result_file
    // std::get<18>(alg_data); // pattern_enumeration_indices
    // std::get<19>(alg_data); // visitor_set  
    auto superstep = std::get<20>(alg_data); // superstep_var
    // std::get<21>(alg_data); // vertex_sequence_number
    // std::get<22>(alg_data); // pattern_aggregation_steps
      
    auto source_index_pattern_indices = 0;

    if (vertex.is_delegate() && (g.master(vertex) != mpi_rank)) {  
      // Important : it should never get here
      //std::cerr << "Error: Controller forwarded visitor to delegates." 
      //  << std::endl;      
      // Important : in the synchronous version it will get here 
//      std::cerr << "Superstep " << superstep 
//        << " Vertex " << g.locator_to_label(vertex) << " "
//        << " is a delegate, ignoring visitor ... " 
//        << std::endl; // Test
      // TODO: update init_visitor_traversal_collection so it ignores delegates?
      return false;
    }
   
//    if (superstep > 0) { // Test
//      std::cerr << "Superstep " << superstep
//        << " Vertex " << g.locator_to_label(vertex) << " "
//        << " visiting ... "
//        << std::endl; // Test 
//    } // Test

    if (ack_success) {
      // token source receiving ack_success 
      auto find_token_source = std::get<6>(alg_data)
        .find(g.locator_to_label(vertex)); // token_source_map
      if (find_token_source == std::get<6>(alg_data).end()) {
        std::cerr << "Error: did not find the expected item in the map." 
          << std::endl;
        //return true; // false ?           
        return false;
      }
      find_token_source->second = 1; //true;   
      std::get<9>(alg_data) = 1; // true; // pattern_found         
//++      return true; // Important : must return true to handle delegates 
      return false;
    } 

    // TODO: remove if delegates do not init/foward tokens 
    // (also, see '++' cpmments 
    // if vertex is a delegate  
    if (!is_init_step && vertex.is_delegate() && 
      (g.master(vertex) != mpi_rank)) { 
      // verify if this vertex has already forwarded the token      
      if (enable_vertex_token_source_cache) {        
      if (!is_init_step && max_itr_count > itr_count) {
        auto find_token_source_forwarded = std::get<12>(alg_data)[vertex]
          .find(g.locator_to_label(target_vertex));
        if (find_token_source_forwarded != 
          std::get<12>(alg_data)[vertex].end()) {
          return false; 
        } 		
      }
      } // enable_vertex_token_source_cache
    }
   
    // TODO: remove if delegates do not init/foward tokens  
    // TODO: verify if this vertex is the 'active' vertex map
    //auto find_vertex = std::get<5>(alg_data).find(g.locator_to_label(vertex));
    //if (find_vertex == std::get<5>(alg_data).end()) {
    //  return false;
    //}    

    if (!do_pass_token && is_init_step && itr_count == 0 && superstep == 0) {
      // init tokens from the source vertex
       
      if (vertex_data != pattern[0]) { 
        return false;  
      }

      // pattern_selected_vertices and vertex_token_source_set 
      if (std::get<16>(alg_data) && std::get<12>(alg_data)[vertex].empty()) { 
        return false;          
      }  

      BitSet vertex_template_vertices(std::get<14>(alg_data)[vertex]); // template_vertices
 
      if (vertex_template_vertices.none() || !vertex_template_vertices.test(pattern_indices[0])) {
        return false;  
      }

      // nonunique metadata
      if (path_checking_filter) { 
        // initiate token only from vertices with multiple template vertex matches 
        if (!expect_target_vertex)  { // path checking
          //if (!(vertex_template_vertices.test(pattern_indices[0]) 
          //  && vertex_template_vertices.test(pattern_indices[pattern_indices.size() -1]) ) ) {
          //  return false;
          //} 
          if (!vertex_template_vertices.test(pattern_indices[0])
            && vertex_template_vertices.count() == 1) {
            return false;
          }
        } 
      } 
      // nonunique metadata    
 
      // TODO: token passing constraints go here
     
      // local, controller and delegates are added to the token_source_map 
      if (!std::get<16>(alg_data)) { // pattern_selected_vertices // the map was populate in the previous iteration    
 
        auto find_token_source = std::get<6>(alg_data).find(g.locator_to_label(vertex)); // token_source_map
        if (find_token_source == std::get<6>(alg_data).end()) {
          auto insert_status = std::get<6>(alg_data).insert({g.locator_to_label(vertex), false});
          if(!insert_status.second) {
            std::cerr << "Error: failed to add an element to the map." << std::endl;
            return false;   
          } 
        }

      }  	
//      }

      // initiate token from the source vertex
//--      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
//--        vertex_locator neighbour = eitr.target();      
      for (auto& item : std::get<15>(alg_data)[vertex]) { // vertex_active_edges_map
        vertex_locator neighbour = g.label_to_locator(item.first);
 
        if (std::get<16>(alg_data)) { // pattern_selected_vertices
          // vertex_token_source_set // the set was populate in the previous iteration 
          for (auto v : std::get<12>(alg_data)[vertex]) {
            tppm_visitor_tds new_visitor(neighbour, vertex, g.label_to_locator(v), 
              g.locator_to_label(neighbour), 
              v, // token_label 
              0, // sequence_number
              visited_vertices, 0, //pattern_cycle_length, 
              //source_index_pattern_indices, 
              pattern_indices[source_index_pattern_indices],
              //pattern_valid_cycle, 
              true, false);
            vis_queue->queue_visitor(new_visitor);
          } // for       
        } else { 
          tppm_visitor_tds new_visitor(neighbour, vertex, vertex, 
            g.locator_to_label(neighbour), 
            g.locator_to_label(vertex), // token_label
            std::get<21>(alg_data)[vertex], //0, // sequence_number 
            visited_vertices, 0, //pattern_cycle_length, 
            //source_index_pattern_indices, 
            pattern_indices[source_index_pattern_indices], 
            // pattern_valid_cycle, 
            true, false);        
          vis_queue->queue_visitor(new_visitor);          
        } 
      }
//++      return true;
      return false; // controller is not forwarding visitor to delegates
 
    } else if (!is_init_step) { // else if     

      // if (superstep > 0) {
      //   itr_count = superstep - 1;
      //   or
      //   new_itr_count = superstep;  
      // }  

      bool do_forward_token = false;
      auto new_itr_count = itr_count + 1; // Important : itr_count is the itr_count of the parent 
      auto next_pattern_index = source_index_pattern_indices + new_itr_count; // expected next pattern_index  
      auto vertex_pattern_index = 0; //find_vertex->second.vertex_pattern_index;
      //auto vertex_pattern_index = pattern_indices[source_index_pattern_indices + new_itr_count]; 
      //// TODO: read from the map
      // vertex_data == pattern[next_pattern_index] 
      // and vertex_pattern_index == pattern_indices[next_pattern_index] must hold  
      // otherwise return false
      if (vertex_data != pattern[next_pattern_index]) {
        return false;
      } 

      bool match_found = false;

      BitSet vertex_template_vertices(std::get<14>(alg_data)[vertex]); // template_vertices

      if (vertex_template_vertices.none()) {
        return false;
      }

      // verify template vertex match  
      // TODO: replace the loop by vertex_template_vertices.test(pattern_indices[next_pattern_index])
      for (size_t i = 0; i < vertex_template_vertices.size(); i++) { 
        if (vertex_template_vertices.test(i)) {
          if (i == pattern_indices[next_pattern_index]) {
            match_found = true;
            vertex_pattern_index = i;
            break; 
          }
        }
      }

      if (!match_found) {
        return false;
      }
  
      if (max_itr_count > itr_count) {
        // verify if received from a valid parent       
        // are vertex_data and vertex_pattern_index valid
//        if (vertex_data == pattern[pattern_indices[next_pattern_index]] && 
	// TODO: Important: for token passing only        
        if (vertex_data == pattern[next_pattern_index] &&
          vertex_pattern_index == pattern_indices[next_pattern_index]) {
          // verify if received from a valid parent
//          for (auto e = pattern_graph.vertices[vertex_pattern_index]; 
//            e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {  
//            if (pattern_graph.edges[e] == parent_pattern_index) {
//              do_forward_token = true; 
//              break; 
//            }  
//          } // for      
          if (parent_pattern_index == pattern_indices[next_pattern_index - 1]) { 
            do_forward_token = true; 
          }   
        } // if    
       
        // TDS 
        //std::get<18>(alg_data); // pattern_enumeration_indices
        if (std::get<18>(alg_data)[new_itr_count] == new_itr_count) {
          // duplicate is not expected in visited_vertices  
          for (size_t i = 0; i < new_itr_count; i++) { // Important : must be < not <=
            if (g.locator_to_label(visited_vertices[i]) == g.locator_to_label(vertex)) {
              //break;
              return false;
            }
          }  
        } else if (std::get<18>(alg_data)[new_itr_count] < new_itr_count) {
          // vertex is expected in visited_vertices
          if (g.locator_to_label(visited_vertices[std::get<18>(alg_data)[new_itr_count]]) 
            != g.locator_to_label(vertex)) {
            return false;
          }
        } else {
          std::cerr << "Error: invalid value." << std::endl;
          return false;
        }         
 
      } else if (max_itr_count == itr_count) {
        // are vertex_data and vertex_pattern_index valid
        //bool match_found = false;
        match_found = false;
 
        // verify if received from a valid parent 
//        if (vertex_data == pattern[pattern_indices[next_pattern_index]] &&
        // TODO: Important: for token passing only        
        if (vertex_data == pattern[next_pattern_index] &&        
          vertex_pattern_index == pattern_indices[next_pattern_index]) { 
          // verify if received from a valid parent
//          for (auto e = pattern_graph.vertices[vertex_pattern_index];
//            e < pattern_graph.vertices[vertex_pattern_index + 1]; e      return false;) {
//            if (pattern_graph.edges[e] == parent_pattern_index) {
//              match_found = true;
//              break;
//            }
//          } // for
          if (parent_pattern_index == pattern_indices[next_pattern_index - 1]) {
            match_found = true; 
          }      
        } // if
        
        if (match_found && !expect_target_vertex) {
          // path checking 
          if (g.locator_to_label(vertex) == g.locator_to_label(target_vertex)) {
           return false; //true; 
                      
          } else {           
            std::get<21>(alg_data)[vertex]++; // update the sequence number    
            // valid path found, send ack to the token source so it could update its state
            tppm_visitor_tds new_visitor(target_vertex, vertex, target_vertex, 
              g.locator_to_label(target_vertex), 
              g.locator_to_label(vertex), //g.locator_to_label(target_vertex), // token_label
              std::get<21>(alg_data)[vertex], //sequence_number, // sequence_number 
              visited_vertices, itr_count, //max_itr_count,
              //source_index_pattern_indices, 
              0, //expect_target_vertex, 
              false, false, true);
            vis_queue->queue_visitor(new_visitor); 

            // write to the paths_result_file
            std::get<17>(alg_data) << "[" << mpi_rank << "], ";
            for (size_t i = 0; i <= new_itr_count ; i++) {
              std::get<17>(alg_data) << g.locator_to_label(visited_vertices[i]) << ", ";
            }
            std::get<17>(alg_data) << "[" << g.locator_to_label(vertex) << "]\n";

            return false; //true;            
          }  
        } else if (g.locator_to_label(vertex) == g.locator_to_label(target_vertex)
          && (g.locator_to_label(vertex) == g.locator_to_label(visited_vertices[0]))
          && match_found && expect_target_vertex) {
          // cycle checking        

          auto find_token_source = std::get<6>(alg_data).find(g.locator_to_label(vertex)); // token_source_map
          if (find_token_source == std::get<6>(alg_data).end()) {
            std::cerr << "Error: did not find the expected item in the map." << std::endl;
            //return true; // false ?           
            return false;
          }

          find_token_source->second = 1; //true;   
 
          std::get<9>(alg_data) = 1; // true; // pattern_found
  
          // write to the paths_result_file  
          std::get<17>(alg_data) << "[" << mpi_rank << "], "; 
          for (size_t i = 0; i <= new_itr_count ; i++) { 
            std::get<17>(alg_data) << g.locator_to_label(visited_vertices[i]) << ", ";  
          }
          std::get<17>(alg_data) << "[" << g.locator_to_label(vertex) << "]\n";          

//++          return true; // Important : must return true to handle delegates
	  return false;
        } else {
          std::cerr << "Error: wrong code branch." << std::endl;
          return false;  
        }   
      } else {
        std::cerr << "Error: wrong code branch." << std::endl;    
        return false;
      }   

      if (!do_forward_token) {
        return false;
      }       
      
      // OK to forward the token

      // if vertex is a delegate // this is important
/*++      if (vertex.is_delegate() && (g.master(vertex) != mpi_rank)) {
        // forwarded a token from a source, now update vertex_token_source_set 
        auto find_token_source_forwarded = std::get<12>(alg_data)[vertex].find(g.locator_to_label(target_vertex));
        if (find_token_source_forwarded == std::get<12>(alg_data)[vertex].end()) {
          auto insert_status = std::get<12>(alg_data)[vertex].insert(g.locator_to_label(target_vertex));		
	  if(!insert_status.second) {
            std::cerr << "Error: failed to add an element to the set." << std::endl;
            return false;
          }
          //std::cout << g.locator_to_label(vertex) << " adding " << g.locator_to_label(target_vertex) 
          //  << " to the vertex set" << std::endl; // Test 
        } else {
          std::cerr << "Error: unexpected item in the set." << std::endl;
	  return false;
        }        
      } // if vertex is a delegate */

      for (auto item : std::get<15>(alg_data)[vertex]) { // vertex_active_edges_map
        vertex_locator neighbour = g.label_to_locator(item.first);
        
 	//std::cout << g.locator_to_label(vertex) << " "
        // << " neighbour count : " << std::get<15>(alg_data)[vertex].size() << " "   
 	// << (new_itr_count + 1) << " " << std::get<18>(alg_data)[new_itr_count + 1] << " "
        // << g.locator_to_label(visited_vertices[std::get<18>(alg_data)[new_itr_count + 1]]) << " "
 	// << g.locator_to_label(neighbour) << std::endl; // Test

        // penultimate hop, only forward to the "valid" target_vertex 
        if (max_itr_count == new_itr_count) { // Important
          if (expect_target_vertex && (g.locator_to_label(target_vertex) != g.locator_to_label(neighbour))) {
            // cycle checking
            continue;
          } else if (!expect_target_vertex && (g.locator_to_label(target_vertex) == g.locator_to_label(neighbour))) {
            // path checking
            continue;
          } else if (!expect_target_vertex && (g.locator_to_label(target_vertex) != g.locator_to_label(neighbour))) {
            // path checking
            match_found = false;
             
            //std::get<18>(alg_data); // pattern_enumeration_indices
            if (std::get<18>(alg_data)[new_itr_count + 1] == new_itr_count + 1) {
              // duplicate is not expected in visited_vertices  
              for (size_t i = 0; i <= new_itr_count; i++) {
                if (g.locator_to_label(visited_vertices[i]) == g.locator_to_label(neighbour)) {
                  match_found = true;
                  break;
                }
              }
            } else if (std::get<18>(alg_data)[new_itr_count + 1] < new_itr_count + 1) {
              // neighbour is expected in visited_vertices
              if (g.locator_to_label(visited_vertices[std::get<18>(alg_data)[new_itr_count + 1]])
                != g.locator_to_label(neighbour)) {
                match_found = true;
              }
            } else {
              std::cerr << "Error: invalid value C." << std::endl;
              //return false;
              match_found = true;
            }

            if (match_found) {
              continue;
            }

          }
        } else if (max_itr_count > itr_count) {
          match_found = false;            
           
          //std::get<18>(alg_data); // pattern_enumeration_indices
          if (std::get<18>(alg_data)[new_itr_count + 1] == new_itr_count + 1) {
            // duplicate is not expected in visited_vertices  
            for (size_t i = 0; i <= new_itr_count; i++) {
              if (g.locator_to_label(visited_vertices[i]) == g.locator_to_label(neighbour)) {                
                match_found = true; 
                break;
              }
            }  
          } else if (std::get<18>(alg_data)[new_itr_count + 1] < new_itr_count + 1) {
            // neighbour is expected in visited_vertices
            if (g.locator_to_label(visited_vertices[std::get<18>(alg_data)[new_itr_count + 1]]) 
              != g.locator_to_label(neighbour)) {
              match_found = true;              
            }  
          } else {
            std::cerr << "Error: invalid value." << std::endl;
            //return false;
            match_found = true;
          }       
           
          if (match_found) {
            continue;
          }      
             
        } else {
          std::cerr << "Error: wrong code branch." << std::endl;
          return false;
        }  

//        std::cerr << "Superstep " << superstep << " "  
//          << g.locator_to_label(vertex) << " vertex_pattern_index " 
//          << vertex_pattern_index << " " << new_itr_count << " forwarding to " 
//          << g.locator_to_label(neighbour) << std::endl; // Test 

        auto new_token_label = target_vertex_label;
        auto new_sequence_number = sequence_number;
        //auto temp_enumeration = false; // Test
        auto enable_new_token_label = std::get<22>(alg_data)[next_pattern_index]; 

        auto has_unique_template_vertex = 
          (vertex_template_vertices.count() == 1 ? true : false);

        /*if (enable_new_token_label && !expect_target_vertex && 
          !has_unique_template_vertex && path_checking_filter) {
          new_token_label = g.locator_to_label(vertex),
          new_sequence_number = std::get<21>(alg_data)[vertex]++; 
        }  

        if (enable_new_token_label && !path_checking_filter) {
          new_token_label = g.locator_to_label(vertex),
          new_sequence_number = std::get<21>(alg_data)[vertex]++;          
        }*/

        if (enable_new_token_label) {
          new_token_label = g.locator_to_label(vertex),
          new_sequence_number = std::get<21>(alg_data)[vertex]++;
        }              
        
        tppm_visitor_tds new_visitor(neighbour, vertex, target_vertex, 
          g.locator_to_label(neighbour), 
          new_token_label, //g.locator_to_label(vertex), //g.locator_to_label(target_vertex), // token_label
          new_sequence_number, //std::get<21>(alg_data)[vertex], //sequence_number, // sequence_number 
          visited_vertices, new_itr_count, //max_itr_count, 
          //source_index_pattern_indices, 
          vertex_pattern_index //, 
          //expect_target_vertex
          ); 
          // vertex_pattern_index = parent_pattern_index for the neighbours 
        vis_queue->queue_visitor(new_visitor);
       
        //if (temp_enumeration) {
        //  std::get<21>(alg_data)[vertex]++; // update the sequence number
        //} 
      }
	     
//++      return true;
      return false; // controller not forwarding visitor to delegates

      // else if
    } else {
      return false;
    }
    return false;		 
  }

  inline bool operator==(const tppm_visitor_tds& visitor) const {     
    return (this->vertex_label == visitor.vertex_label) &&
      (this->target_vertex_label == visitor.target_vertex_label) && 
      (this->sequence_number == visitor.sequence_number);  
  }

  friend inline bool operator>(const tppm_visitor_tds& v1, const tppm_visitor_tds& v2) {
    return false;
    /*if (v1.itr_count > v2.itr_count) {
      return true;
    } else if (v1.itr_count < v2.itr_count) {
      return false; 
    }
    if (v1.vertex == v2.vertex) {
      return false;
    }
    return !(v1.vertex < v2.vertex);*/
    
    /*if (v1.itr_count <= v2.itr_count) {
      return true;
    } else {
      return false;
    }*/ 
  }

  //friend inline bool operator<(const tppm_visitor& v1, const tppm_visitor& v2) {
    //re://www.youtube.com/watch?v=S5fn1DfqPfA&list=RDEMQETG-l6Y1rVZh5nSwVP8Jg&index=13turn false;
    //if (v1.itr_count < v2.itr_count) {
    //  return true;
    //} else if (v1.itr_count > v2.itr_count) {
    //  return false;
    //}
    //if (v1.vertex == v2.vertex) {
    //  return false;
    //}
    //return !(v1.vertex < v2.vertex);    
  //}  

  vertex_locator vertex;
  vertex_locator parent;  
  vertex_locator target_vertex; // for a cycle, this is also the originating vertex
  //size_t source_index_pattern_indices; // index of the token source in the pattern_indices container 
  size_t parent_pattern_index; // TODO: change to the same type as in the pattern_graph // TODO: remove
  std::array<vertex_locator, 13> visited_vertices; // I think we will have to replace vertex_locator with Vertex

  Vertex vertex_label;
  Vertex target_vertex_label; // TODO: this label will be updated at the points where no caching happens // token_label?  
  uint64_t sequence_number;

  size_t itr_count; // TODO: change type // Important : itr_count of the parent // TODO: remove 
  //size_t max_itr_count; // equal to diameter - 1 of the pattern as itr_count is initialized to 0 // TODO: change type

  //bool expect_target_vertex;
  bool do_pass_token; // TODO: remove
  bool is_init_step;
  bool ack_success; 
  
  // TODO: replace all of the above by msg_type
  uint8_t msg_type; // 0 - init, 1 - foward, 2 - ack-success, 3 - ack-failure 
};

template <typename TGraph, typename Vertex, typename Edge, typename VertexData, 
  typename EdgeData, typename VertexMetadata, typename EdgeMetadata, 
  typename VertexActive, typename VertexUint8MapCollection, 
  typename TemplateVertex, typename VertexStateMap, typename PatternGraph, 
  typename PatternUtilities, typename VertexUint8Map, 
  typename VertexSetCollection, 
  template<typename> class DelegateGraphVertexDataSTDAllocator,
  typename Boolean, typename BitSet>      

void token_passing_pattern_matching(TGraph* g, VertexMetadata& vertex_metadata,
  VertexActive& vertex_active, 
  VertexUint8MapCollection& vertex_active_edges_map, 
  TemplateVertex& template_vertices, VertexStateMap& vertex_state_map,
  PatternGraph& pattern_graph, PatternUtilities& pattern_utilities, size_t pl, 
  VertexUint8Map& token_source_map, 
  VertexSetCollection& vertex_token_source_set,
  std::vector<uint8_t>::reference pattern_found, 
  std::ofstream& paths_result_file, uint64_t& message_count) {

  int mpi_rank = havoqgt_env()->world_comm().rank();
  size_t superstep_var = 0;
  size_t& superstep_ref = superstep_var;  

  if (mpi_rank == 0) {
    std::cout << "Token Passing ... " << std::endl;
  }

  // TODO: temporary patch
  auto pattern = std::get<0>(pattern_utilities.input_patterns[pl]);
  auto pattern_indices = std::get<1>(pattern_utilities.input_patterns[pl]);
  auto pattern_cycle_length = std::get<2>(pattern_utilities.input_patterns[pl]); // uint
  auto pattern_valid_cycle = std::get<3>(pattern_utilities.input_patterns[pl]); // boolean
  auto pattern_interleave_label_propagation = std::get<4>(pattern_utilities.input_patterns[pl]); // boolean
//--  auto pattern_seleted_edges_tp = std::get<5>(ptrn_util_two.input_patterns[pl]); // boolean 
  auto pattern_selected_vertices = std::get<5>(pattern_utilities.input_patterns[pl]); // boolean
  auto pattern_enumeration_indices = pattern_utilities.enumeration_patterns[pl];
  auto pattern_aggregation_steps = pattern_utilities.aggregation_steps[pl];
  // TODO: temporary patch
  
  // TODO: remove
  uint8_t vertex_rank = 0; // dummy 
  uint8_t edge_metadata = 0; // dummy 
  // TODO: remove   
   
   
  uint64_t max_itr_count = std::get<2>(pattern_utilities.input_patterns[pl]);
  //uint64_t itr_count = 0; // superstep_var
  bool expect_target_vertex = std::get<3>(pattern_utilities.input_patterns[pl]);
  //Vertex parent_pattern_vertex = ; // a.k.a parent_pattern_index;
  auto source_index_pattern_indices = 0;    

  bool enable_filtered_token_relay = false; 

  typedef tppm_visitor_tds<TGraph, Vertex, BitSet> visitor_type;

  typedef std::vector<visitor_type> VectorVisitor;
  typedef DelegateGraphVertexDataSTDAllocator<VectorVisitor> VertexVisitorCollection;
  VertexVisitorCollection vertex_visitors(*g);

  typedef DelegateGraphVertexDataSTDAllocator<uint64_t> VertexUint64Collection;
  VertexUint64Collection vertex_sequence_number(*g);

  struct VisitorCompare {
    public:
      bool operator()(const visitor_type& v1, const visitor_type& v2) const {
        auto return_value = (v1.vertex_label == v2.vertex_label) && 
          (v1.target_vertex_label == v2.target_vertex_label) && 
          (v1.sequence_number == v2.sequence_number);   
        return return_value; 
      }			
  }; 

  struct VisitorHash {
    public:
      size_t operator()(const visitor_type& visitor) const {        
        return visitor.vertex_label;
      }
  };

  typedef std::unordered_set<visitor_type, VisitorHash, VisitorCompare> VisitorSet;
  VisitorSet* visitor_set_receive = new VisitorSet();
  VisitorSet* visitor_set_send = new VisitorSet();  
  visitor_set_receive->clear();
  visitor_set_send->clear(); 

  auto alg_data = std::forward_as_tuple(vertex_metadata, pattern, pattern_indices, vertex_rank, 
    pattern_graph, vertex_state_map, token_source_map, pattern_cycle_length, pattern_valid_cycle, pattern_found, 
    edge_metadata, g, vertex_token_source_set, vertex_active, template_vertices, vertex_active_edges_map, 
    pattern_selected_vertices, paths_result_file, pattern_enumeration_indices, visitor_set_receive, 
    superstep_var, vertex_sequence_number, pattern_aggregation_steps);

  auto vq = create_visitor_queue<visitor_type, /*havoqgt::detail::visitor_priority_queue*/tppm_queue>(g, alg_data);
 
  for (auto superstep = 0; superstep <= max_itr_count /*<=pattern_cycle_length*/; superstep++) {
    superstep_ref = superstep;
     
    if (mpi_rank == 0) {
      std::cout << "Token Passing | Superstep #" << superstep;
    }       

    double time_start = MPI_Wtime();
    
    if (superstep == 0) {
      if (mpi_rank == 0) {
        std::cout << " | Initiating ..."; // Test    
      } 
      vq.init_visitor_traversal_new(); 
      //std::swap(visitor_set_receive, visitor_set_send);
      //visitor_set_receive->clear();     
    } else {       
      if (mpi_rank == 0) {
        std::cout << " | Forwarding ..."; // Test  
      } 
      //std::swap(visitor_set_receive, visitor_set_send);
      //visitor_set_receive->clear();
//      if (visitor_set_send->size() > 0) { 
//        std::cout << "\nMPI Rank : " << mpi_rank << " Superstep " << superstep 
//          << " visitor set send size " << visitor_set_send->size() << std::endl; // Test
//      }

      vq.init_visitor_traversal_collection(*visitor_set_send);
       
      // delegates must ignore the visitors 
      //std::swap(visitor_set_receive, visitor_set_send); 
      //visitor_set_receive->clear();  
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout <<  " | Synchronizing ... "; 
    }
  
//    if (visitor_set_receive->size() > 0) {
//    std::cout << "\nMPI Rank : " << mpi_rank << " Superstep " << superstep
//      << " visitor set receive size " << visitor_set_receive->size() << std::endl; // Test
//    }

    size_t global_visitor_set_size =
      havoqgt::mpi::mpi_all_reduce(visitor_set_receive->size(),  
      std::plus<size_t>(), MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << " | Global Visitor Count : "
        << global_visitor_set_size;
    }
   
    std::swap(visitor_set_receive, visitor_set_send);
    visitor_set_receive->clear();
    vertex_sequence_number.reset(0);

    MPI_Barrier(MPI_COMM_WORLD); // TODO: ?

    double time_end = MPI_Wtime();
    if (mpi_rank == 0) {
      std::cout << " | Time : " << time_end - time_start << std::endl;  
    }      

  } // for

  visitor_set_receive->clear();
  visitor_set_send->clear();

}

}} //end namespace havoqgt::mpi

#endif //HAVOQGT_TOKEN_PASSING_PATTERN_MATCHING_NONUNIQUE_ITERATIVE_TDS_1_HPP_INCLUDED
