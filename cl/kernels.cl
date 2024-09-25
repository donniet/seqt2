
// index of first element in sorted that is >= target
long lower_bound(long target, global long * sorted, long length) {
    long low = 0;
    long high = length - 1;

    while(low <= high) {
        long i = (low + high) / 2;

        if(sorted[i] >= target) 
            high = i - 1;
        else
            low = i + 1;
    } 

    return low;
}

// index of first element in sorted that is >= target

bool greater_equal(long2 x, long2 y) {
    if(x.s0 > y.s0)
        return true;
    if(x.s0 < y.s0)
        return false;
    if(x.s1 >= y.s1)
        return true;

    return false;
}

bool equal(long2 x, long2 y) {
    return x.s0 == y.s0 && x.s1 == y.s1;
}

long lower_bound2(long2 target, global long2 * sorted, long length) {
    long low = 0;
    long high = length - 1;

    while(low <= high) {
        long i = (low + high) / 2;

        if(greater_equal(sorted[i], target))
            high = i - 1;
        else
            low = i + 1;
    } 

    return low;
}

// index of first element in sorted that is > target
long upper_bound(long target, global long * sorted, long length) {
    long low = 0;
    long high = length;

    while(low < high) {
        long i = (low + high) / 2;

        if(target >= sorted[i]) 
            low = i + 1;
        else
            high = i;
    } 

    return low;
}

kernel void is_sequence_significant(
    global long2 * seq,
    long total_sequences,
    global float * significance,
    global long * counts,
    long total,
    float sigma,
    long min_count,
    global long * output
) {
    const long gid = get_global_id(0);
    if(gid >= total_sequences)
        return;

    if(significance[seq[gid].s0] > sigma && significance[seq[gid].s1] > sigma && 
        counts[seq[gid].s0] > min_count && counts[seq[gid].s1] > min_count)
        output[gid] = 1;
    else
        output[gid] = 0;
}

kernel void find_nexts(
    global long * tracked, 
    long total,
    global long * sorted_current_lengths,
    long current_total,
    global long2 * nexts)
{
    const long gid = get_global_id(0);
    if(gid >= total)
        return;

    // is the previous one behind by the length of the current ones?
    long search_length = -tracked[gid];

    // do a binary search through the sorted currents for search_length
    nexts[gid].s0 = lower_bound(search_length, sorted_current_lengths, current_total);
    nexts[gid].s1 = upper_bound(search_length, sorted_current_lengths, current_total);
    // nexts_begin[gid] = 0;
    // nexts_end[gid] = 1;
}

kernel void gather_seqs(
    global const long * index,
    long index_size,
    global const long2 * seqs,
    global const long * reverse_index,
    global long2 * output
) {
    const long gid = get_global_id(0);
    if(gid >= index_size)
        return;

    long2 seq = seqs[index[gid]];

    seq.s0 = reverse_index[seq.s0];
    seq.s1 = reverse_index[seq.s1];

    output[gid] = seq;
}

// TODO: this could be way more efficient
kernel void depends_on_flagged(
    global const long2 * seqs,
    long seqs_size,
    global long * flagged,
    global long * output
) {
    const long gid = get_global_id(0);
    if(gid >= seqs_size)
        return;
    
    if(flagged[gid] || flagged[seqs[gid].s0] || flagged[seqs[gid].s1])
        output[gid] = 1;
    else
        output[gid] = 0;
}

kernel void depends_on_sorted_list(
    global const long2 * seqs,
    global const long * sorted,
    long sorted_total,
    global const long * begin,
    long total,
    global long * output
) {
    const long gid = get_global_id(0);
    if(gid >= total)
        return;

    long check = seqs[begin[gid]].s0;

    long low = lower_bound(check, sorted, sorted_total);
    long upp = upper_bound(check, &sorted[low], sorted_total - low);

    if(low < upp) {
        output[gid] = 1;
        return;
    }

    check = seqs[begin[gid]].s1;

    low = lower_bound(check, sorted, sorted_total);
    upp = upper_bound(check, &sorted[low], sorted_total - low);

    if(low < upp) {
        output[gid] = 1;
        return;
    }

    output[gid] = 0;
}

kernel void calculate_stats(
    global long2 * seqs,
    global long2 * initial_seq_counts,
    global long * counts,
    global long * initial_characters_read,
    global long * lengths,
    long characters_read,
    global float * stddev_counts,
    global float * expected_counts,
    global float * significance,
    long total
) {
    const long gid = get_global_id(0);
    if(gid >= total)
        return;

    // how many a's and b's have we seen since we started tracking ab?
    long a = counts[seqs[gid].s0] - initial_seq_counts[gid].s0;
    long b = counts[seqs[gid].s1] - initial_seq_counts[gid].s1;

    if(seqs[gid].s0 == 0 || seqs[gid].s1 == 0 || a == 0 || b == 0) {
        // this is either an atom which is as significant as it's count, or it's a newly initialized sequence
        stddev_counts[gid] = 1;
        expected_counts[gid] = 0;
        significance[gid] = (float)counts[gid];
        return;
    } 

    long a_len = lengths[seqs[gid].s0];
    long b_len = lengths[seqs[gid].s1];
    long min_initial = min(initial_characters_read[seqs[gid].s0], initial_characters_read[seqs[gid].s1]);
    float characters_since = (float)(characters_read - min_initial);

    // characters since a was first spotted minus it's own characters
    // we subtract one from the overall count because we start initial_characters after the first occurance
    float characters_since_a = (float)(characters_read - initial_characters_read[seqs[gid].s0] - a_len * (a - 1));

    // characters since b was first completed minus it's own characters
    float characters_since_b = (float)(characters_read - initial_characters_read[seqs[gid].s1] - b_len * (b - 1));

    // probability that at a given location what follows is a
    float pa = characters_since_a == 0 ? 1e10 : (float)a / characters_since_a; // As/chars
    // probability that at the given location what is behind is b
    float pb = characters_since_b == 0 ? 1e10 : (float)b / characters_since_b; // Bs/chars
    // probability that we are right between a sequence of ab
    float pab = pa * pb; // As*Bs/chars^2



    // probability that any particular location is between a and b
    expected_counts[gid] = (characters_since - b_len) * pab; // As * Bs / chars

    if(characters_since - b_len == 0) {
        stddev_counts[gid] = 1;
        significance[gid] = (float)counts[gid];
        return;
    }

    float variation = (characters_since - b_len) * pab * (1. - pab); // As^2 * Bs^2 / chars^3

    stddev_counts[gid] = sqrt(variation);

    significance[gid] = ((float)counts[gid] - expected_counts[gid]) / stddev_counts[gid];
}

kernel void collect_finds(
    global long2 * nexts,
    global long * scratch,
    long total,
    global long * current,
    global long2 * found,
    long found_count
) {
    const long gid = get_global_id(0);
    if(gid >= found_count)
        return;

    // gid is the index of this find from the sum of scratch

    // search through the sorted scratch vector to get the index of prev

    // the upper_bound of gid will give us the subsequent 
    // because gid ranges continuously from [0,found_count)
    // and scratch is an inclusive_scan with each step being
    // the length of the next finds
    long prev = upper_bound(gid, scratch, total);
    
    // prev == 0 means that index 0 (our begin token)
    // has at least one next (scratch[0] > 0 where
    // scratch is a scan of the lengths of nexts)
    // in this case our "first_find" should be 0
    long first_find = 0;
    if(prev > 0) 
        first_find = scratch[prev-1];
    
    // our next is the offset of where our nexts begin
    // by our gid from the first find
    long next = current[nexts[prev].s0 + gid - first_find];

    // output our results
    found[gid].s0 = prev;
    found[gid].s1 = next;
}

kernel void scatter_value(
    global long * indices,
    long total_indices,
    long value,
    global long * output)
{
    const long gid = get_global_id(0);
    if(gid >= total_indices)
        return;

    output[indices[gid]] = value;
}

kernel void mark_exists(
    global long2 * seqs,
    long total,
    global long2 * found,
    global long * scratch,
    long found_count
) {
    const long gid = get_global_id(0);
    if(gid >= total)
        return;

    // we are looking at all our existing nodes and marking
    // scratch if this one exists in found
    long2 value = seqs[gid];

    long low = lower_bound2(value, found, found_count);
    if(equal(found[low], value)) {
        scratch[low] = gid;
    }
}

kernel void initialize_newly_found_sequences(
    global long * new_indices,
    long new_total,
    global long2 * found,
    global long * lengths,
    global long2 * seqs,
    global long2 * initial_seq_counts,
    global long * tracked,
    global long * counts,
    long tracked_value,
    long previous_total
) {
    const long gid = get_global_id(0);
    if(gid >= new_total) 
        return;

    long2 s = found[new_indices[gid]];

    long new_index = previous_total + gid;
    lengths[new_index] = lengths[s.s0] + lengths[s.s1];
    seqs[new_index] = s;
    counts[new_index] = 0;
    initial_seq_counts[new_index].s0 = counts[s.s0];
    initial_seq_counts[new_index].s1 = counts[s.s1];
    tracked[new_index] = tracked_value;
}

kernel void make_pair_constant_second(global long * first, long second, global long2 * output, long total) {
    const long gid = get_global_id(0);
    if(gid >= total)
        return;

    output[gid].s0 = first[gid];
    output[gid].s1 = second;
}

kernel void pack(long total, global long *scanned, global long *output) {
    const long gid = get_global_id(0);
    if (gid >= total) {
        return;
    }

    if (gid == 0) {
        if (scanned[0] != 0) {
            output[0] = 0;
        }
    } else if (scanned[gid] > scanned[gid - 1]) {
        output[scanned[gid - 1]] = gid;
    }
}
