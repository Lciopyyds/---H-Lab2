#include <cstdio>
#include <vector>
#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <iostream>
#include <tuple>
#include <fstream>

// 定义节点类型和位置信息
struct AlignmentNode {
    enum Type { CONTROL, FORWARD, REVERSE } type;
    unsigned int reference_pos;
    unsigned int query_pos;
    
    bool operator==(const AlignmentNode& other) const {
        return type == other.type && 
               reference_pos == other.reference_pos && 
               query_pos == other.query_pos;
    }
};

// 哈希函数特化
namespace std {
template <>
struct hash<AlignmentNode> {
    size_t operator()(const AlignmentNode& node) const {
        return (static_cast<size_t>(node.type) << 32) ^ 
               (node.reference_pos << 16) ^ 
               node.query_pos;
    }
};
}

// 节点距离信息
struct NodeDistance {
    int distance;
    AlignmentNode previous_node;
};

class SequenceAligner {
private:
    std::vector<char> reference_sequence;
    std::vector<char> reverse_complement;
    std::vector<char> query_sequence;
    
    std::unordered_map<AlignmentNode, NodeDistance> distance_map;
    std::unordered_set<AlignmentNode> visited_nodes;
    std::deque<AlignmentNode> processing_queue;
    
    static constexpr int JUMP_WINDOW_SIZE = 1200;
    static constexpr int MIN_KMER_LENGTH = 10;

public:
    SequenceAligner() = default;
    
    void loadSequences(const std::string& query_path, const std::string& reference_path) {
        std::ifstream query_file(query_path);
        std::ifstream ref_file(reference_path);
        
        if (!query_file.is_open() || !ref_file.is_open()) {
            throw std::runtime_error("Failed to open input files.");
        }

        std::string query_str, ref_str;
        std::getline(query_file, query_str);
        std::getline(ref_file, ref_str);

        // 处理参考序列
        for (char nucleotide : ref_str) {
            reference_sequence.push_back(nucleotide);
            reverse_complement.push_back(complement(nucleotide));
        }
        
        // 处理查询序列
        for (char nucleotide : query_str) {
            query_sequence.push_back(nucleotide);
        }
    }
    
    void findOptimalAlignment() {
        initializeSearch();
        performGraphSearch();
        outputAlignmentResults();
    }

private:
    static char complement(char nucleotide) {
        switch (nucleotide) {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'C': return 'G';
            case 'G': return 'C';
            default: return nucleotide;
        }
    }
    
    void initializeSearch() {
        AlignmentNode start_node{AlignmentNode::CONTROL, 0, 0};
        distance_map[start_node] = {0, start_node};
        processing_queue.push_back(start_node);
    }
    
    void performGraphSearch() {
        const unsigned int target_query_pos = query_sequence.size();
        
        while (!processing_queue.empty()) {
            AlignmentNode current_node = processing_queue.front();
            
            if (current_node.query_pos == target_query_pos) {
                std::cout << "Find a path with distance of " << distance_map[current_node].distance << std::endl;
                break;
            }
            
            processing_queue.pop_front();
            
            if (visited_nodes.count(current_node)) continue;
            visited_nodes.insert(current_node);
            
            int current_distance = distance_map[current_node].distance;
            
            switch (current_node.type) {
                case AlignmentNode::CONTROL:
                    processControlNode(current_node, current_distance);
                    break;
                case AlignmentNode::FORWARD:
                    processForwardNode(current_node, current_distance);
                    break;
                case AlignmentNode::REVERSE:
                    processReverseNode(current_node, current_distance);
                    break;
            }
        }
    }
    
    void processControlNode(const AlignmentNode& node, int current_distance) {
        // 计算跳转窗口边界
        unsigned int window_start = (node.query_pos > JUMP_WINDOW_SIZE) ? 
                                    node.query_pos - JUMP_WINDOW_SIZE : 0;
        unsigned int window_end = (node.query_pos + JUMP_WINDOW_SIZE < reference_sequence.size()) ? 
                                  node.query_pos + JUMP_WINDOW_SIZE : reference_sequence.size();
        
        for (unsigned int ref_pos = window_start; ref_pos <= window_end; ref_pos++) {
            // 前向比对节点
            AlignmentNode forward_node{AlignmentNode::FORWARD, ref_pos, node.query_pos};
            updateNodeDistance(forward_node, current_distance + 1, node);
            
            // 反向比对节点
            AlignmentNode reverse_node{AlignmentNode::REVERSE, ref_pos, node.query_pos};
            updateNodeDistance(reverse_node, current_distance + 1, node);
        }
        
        
    }
    
    void processForwardNode(const AlignmentNode& node, int current_distance) {
        // 处理匹配/不匹配情况
        if (node.reference_pos + 1 <= reference_sequence.size() && 
            node.query_pos + 1 <= query_sequence.size()) {
            
            AlignmentNode next_node{AlignmentNode::FORWARD, node.reference_pos + 1, node.query_pos + 1};
            int cost = (reference_sequence[node.reference_pos] == query_sequence[node.query_pos]) ? 0 : 1;
            updateNodeDistance(next_node, current_distance + cost, node, cost == 0);
        }
        
        // 处理参考序列插入
        if (node.reference_pos + 1 <= reference_sequence.size()) {
            AlignmentNode next_node{AlignmentNode::FORWARD, node.reference_pos + 1, node.query_pos};
            updateNodeDistance(next_node, current_distance + 1, node);
        }
        
        // 处理查询序列插入
        if (node.query_pos + 1 <= query_sequence.size()) {
            AlignmentNode next_node{AlignmentNode::FORWARD, node.reference_pos, node.query_pos + 1};
            updateNodeDistance(next_node, current_distance + 1, node);
        }
        
        // 跳转回控制节点
        AlignmentNode control_node{AlignmentNode::CONTROL, 0, node.query_pos};
        updateNodeDistance(control_node, current_distance + 1, node);
    }
    
    void processReverseNode(const AlignmentNode& node, int current_distance) {
        // 处理反向匹配/不匹配情况
        if (node.reference_pos >= 1 && node.query_pos + 1 <= query_sequence.size()) {
            AlignmentNode next_node{AlignmentNode::REVERSE, node.reference_pos - 1, node.query_pos + 1};
            int cost = (reverse_complement[node.reference_pos - 1] == query_sequence[node.query_pos]) ? 0 : 1;
            updateNodeDistance(next_node, current_distance + cost, node, cost == 0);
        }
        
        // 处理反向参考序列插入
        if (node.reference_pos >= 1) {
            AlignmentNode next_node{AlignmentNode::REVERSE, node.reference_pos - 1, node.query_pos};
            updateNodeDistance(next_node, current_distance + 1, node);
        }
        
        // 处理反向查询序列插入
        if (node.query_pos + 1 <= query_sequence.size()) {
            AlignmentNode next_node{AlignmentNode::REVERSE, node.reference_pos, node.query_pos + 1};
            updateNodeDistance(next_node, current_distance + 1, node);
        }
        
        // 跳转回控制节点
        AlignmentNode control_node{AlignmentNode::CONTROL, 0, node.query_pos};
        updateNodeDistance(control_node, current_distance + 1, node);
    }
    
    void updateNodeDistance(const AlignmentNode& node, int new_distance, 
                           const AlignmentNode& previous_node, bool high_priority = false) {
        auto it = distance_map.find(node);
        if (it == distance_map.end() || it->second.distance > new_distance) {
            distance_map[node] = {new_distance, previous_node};
            if (high_priority) {
                processing_queue.push_front(node);
            } else {
                processing_queue.push_back(node);
            }
        }
    }
    
    void outputAlignmentResults() {
        if (processing_queue.empty()) return;
        
        std::vector<std::tuple<int, int, int, int>> aligned_segments;
        AlignmentNode current_node = processing_queue.front();
        int segment_end_ref = current_node.reference_pos;
        int segment_end_query = current_node.query_pos;
        
        // 回溯路径
        while (!(current_node.type == AlignmentNode::CONTROL && 
                 current_node.reference_pos == 0 && 
                 current_node.query_pos == 0)) {
            
            const auto& prev_info = distance_map[current_node];
            
            // 记录比对片段
            if ((current_node.type == AlignmentNode::FORWARD || 
                 current_node.type == AlignmentNode::REVERSE) && 
                current_node.type != prev_info.previous_node.type) {
                
                int query_start = current_node.query_pos;
                int query_end = segment_end_query;
                int ref_start = (current_node.type == AlignmentNode::FORWARD) ? 
                                 current_node.reference_pos : segment_end_ref - 1;
                int ref_end = (current_node.type == AlignmentNode::FORWARD) ? 
                               segment_end_ref : current_node.reference_pos - 1;
                
                if (query_end - query_start + 1 > MIN_KMER_LENGTH) {
                    aligned_segments.emplace_back(query_start, query_end, ref_start, ref_end);
                }
            }
            
            // 更新片段结束位置
            if (current_node.type == AlignmentNode::CONTROL && 
                current_node.type != prev_info.previous_node.type) {
                segment_end_ref = prev_info.previous_node.reference_pos;
                segment_end_query = prev_info.previous_node.query_pos;
            }
            
            current_node = prev_info.previous_node;
        }
        
        // 输出比对结果
        std::cout << "[";
        for (size_t i = aligned_segments.size(); i > 0; --i) {
            const auto& [q_start, q_end, r_start, r_end] = aligned_segments[i - 1];
            std::cout << "(" << q_start << "," << q_end << "," << r_start << "," << r_end << ")";
            if (i > 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
    }
};

int main() {
    try {
        SequenceAligner aligner;
        aligner.loadSequences("query2.txt", "reference2.txt");
        aligner.findOptimalAlignment();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
}