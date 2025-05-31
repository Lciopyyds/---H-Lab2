# 算法设计与分析（H）---Lab2

## 算法伪代码

```c
CONSTANTS:
    FORWARD = "FWD"
    REVERSE = "REV"
    CONTROL = "CTRL"
    MAX_JUMP = 1000
    MIN_SEG_LEN = 15

CLASS SegmentNode:
    kind: String  // FORWARD, REVERSE, or CONTROL
    ref_idx: Int  // Reference index
    qry_idx: Int  // Query index

CLASS PathInfo:
    cost: Int     // Cost to reach this node
    parent: SegmentNode  // Parent node in the path

FUNCTION comp(base):
    RETURN complementary base of input base

CLASS GenomeMatcher:
    query: String
    reference: String
    ref_rev: String  // Reverse complement of reference
    
    memo: Map<SegmentNode, PathInfo>  // Memoization table
    processed: Set<SegmentNode>      // Processed nodes
    queue: Deque<SegmentNode>         // Processing queue (double-ended)
    
    CONSTRUCTOR(query, reference):
        this.query = query
        this.reference = reference
        ref_rev = reversed_complement(reference)
    
    // Main alignment function
    FUNCTION align():
        initialize()
        search_paths()
        print_results()
    
    FUNCTION initialize():
        start = SegmentNode(CONTROL, 0, 0)
        memo[start] = PathInfo(0, start)
        queue.push_back(start)
    
    FUNCTION search_paths():
        end_pos = query.length()
        WHILE not queue.empty():
            current = queue.pop_front()
            
            // Reached end of query sequence
            IF current.qry_idx == end_pos:
                queue.push_front(current)  // Save for backtracking
                BREAK
            
            IF current in processed:
                CONTINUE
            processed.add(current)
            dist = memo[current].cost
            
            SWITCH current.kind:
                CASE CONTROL:
                    expand_control(current, dist)
                CASE FORWARD:
                    expand_forward(current, dist)
                CASE REVERSE:
                    expand_reverse(current, dist)
    
    // Expand CONTROL state
    FUNCTION expand_control(node, dist):
        q = node.qry_idx
        start = max(0, q - MAX_JUMP)
        end = min(reference.length, q + MAX_JUMP)
        FOR r = start TO end:
            try_update(FORWARD, r, q, dist + 1, node)
            try_update(REVERSE, r, q, dist + 1, node)
    
    // Expand FORWARD state
    FUNCTION expand_forward(node, dist):
        r = node.ref_idx
        q = node.qry_idx
        // Match next nucleotides
        IF r < reference.length AND q < query.length:
            penalty = (reference[r] == query[q]) ? 0 : 1
            try_update(FORWARD, r + 1, q + 1, dist + penalty, node, urgent=(penalty==0))
        // Skip reference nucleotide
        IF r < reference.length:
            try_update(FORWARD, r + 1, q, dist + 1, node)
        // Skip query nucleotide
        IF q < query.length:
            try_update(FORWARD, r, q + 1, dist + 1, node)
        // Return to CONTROL state
        try_update(CONTROL, 0, q, dist + 1, node)
    
    // Expand REVERSE state
    FUNCTION expand_reverse(node, dist):
        r = node.ref_idx
        q = node.qry_idx
        // Match next nucleotides
        IF r > 0 AND q < query.length:
            penalty = (ref_rev[r - 1] == query[q]) ? 0 : 1
            try_update(REVERSE, r - 1, q + 1, dist + penalty, node, urgent=(penalty==0))
        // Skip reference nucleotide
        IF r > 0:
            try_update(REVERSE, r - 1, q, dist + 1, node)
        // Skip query nucleotide
        IF q < query.length:
            try_update(REVERSE, r, q + 1, dist + 1, node)
        // Return to CONTROL state
        try_update(CONTROL, 0, q, dist + 1, node)
    
    FUNCTION try_update(kind, ref_idx, qry_idx, new_cost, parent, urgent=false):
        node = SegmentNode(kind, ref_idx, qry_idx)
        // Update if new path is better or first time visiting this node
        IF node not in memo OR memo[node].cost > new_cost:
            memo[node] = PathInfo(new_cost, parent)
            IF urgent:
                queue.push_front(node)  // High priority
            ELSE:
                queue.push_back(node)  // Normal priority
    
    FUNCTION print_results():
        IF queue.empty():
            PRINT "[]"
            RETURN
        
        node = queue.front()
        segments = []
        q_end = node.qry_idx
        r_end = node.ref_idx
        
        // Backtrack to reconstruct path
        WHILE !(node.kind == CONTROL AND node.ref_idx == 0 AND node.qry_idx == 0):
            prev = memo[node].parent
            
            // Found a segment boundary
            IF (node.kind == FORWARD OR node.kind == REVERSE) AND node.kind != prev.kind:
                q_start = node.qry_idx
                // Calculate segment positions
                IF node.kind == FORWARD:
                    r_start = node.ref_idx
                ELSE:
                    r_start = r_end - 1
                // Filter out short segments
                IF q_end - q_start + 1 > MIN_SEG_LEN:
                    ADD (q_start, q_end, r_start, r_end) to segments
            
            // Save control points for segment boundaries
            IF node.kind == CONTROL AND node.kind != prev.kind:
                r_end = prev.ref_idx
                q_end = prev.qry_idx
            
            node = prev  // Move to parent
        
        // Output segments in correct order
        PRINT "["
        FOR i = segments.length - 1 DOWNTO 0:
            PRINT "({segments[i].q_start},{segments[i].q_end},{segments[i].r_start},{segments[i].r_end})"
            IF i > 0: PRINT ", "
        PRINT "]"

MAIN FUNCTION:
    READ query and reference sequences
    CREATE GenomeMatcher(query, reference)
    CALL matcher.align()
```

## 算法复杂度

### 时间复杂度

该算法的时间复杂度为 $O(n\cdot m\cdot MAX_{JUMP})$，其中：

- n = 查询序列长度 (query)
- m = 参考序列长度 (reference)
- MAX_JUMP = 最大跳跃距离 (常量，默认1000)

主要复杂度来源：

1. **状态空间**：每个状态由(qry_idx, ref_idx, kind)三个参数定义，复杂度为 O(n*m*3)
2. **控制状态扩展**：每个控制状态扩展 O(MAX_JUMP) 个子状态
3. **优先队列操作**：每个状态都会被处理一次，处理复杂度为 O(MAX_JUMP + log V)，其中V是总状态数

### 空间复杂度

空间复杂度为 **O(n\*m)**，主要来自于：

- memo数据结构存储所有访问过的状态和路径信息
- processed集合记录已处理状态
- 队列存储待处理状态

在极端情况下（完全匹配），空间复杂度可以优化到 O(n + m)，但实际应用中是 O(n*m)
