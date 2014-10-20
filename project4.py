"""
Project 4 - Computing Global and Local Alignments
"""

#alignment matrix will be a list of lists
#scoring matrix will be a dictionary of dictionaries

def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    Input:  set of characters alphabet, integers diag_score, off_diag_score, \
    and dash_score.
    Output:  dictionary of dictionaries assigning a score for each pairing \
    of characters in the alphabet.
    """
    # add a dash to the alphabet
    alphabet_with_dash = alphabet.union('-')
    
    # create dictionary; score accessed by scoring_matrix[row_char][col_char]
    scoring_matrix = {row_char: {} for row_char in alphabet_with_dash}
    for row_char in alphabet_with_dash:
        for col_char in alphabet_with_dash:
            if row_char == '-' or col_char == '-':
                score = dash_score
            elif row_char == col_char:
                score = diag_score
            else:
                score = off_diag_score
            scoring_matrix[row_char][col_char] = score
    return scoring_matrix

#s_matrix = build_scoring_matrix(set(['a', 't', 'c', 'g']), 10, 4, -6)
#s_matrix = build_scoring_matrix(set(['a', 't']), 10, 4, -6)
# print s_matrix

def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    Input:  two strings seq_x and seq_y, a dictionary of dictionaries \
    scoring_matrix, and a boolean global_flag.
    Output:  global alignment grid represented as a list of lists.
    """
    # initialize empty alignment matrix, len(seq_x) by len(seq_y)
    alignment_matrix = [[0 for dummy_y in range(len(seq_y)+1)] 
                        for dummy_x in range(len(seq_x)+1)]
    
    # compute the score per Q8 in the homework
    # if global flag is True and computed val is negative -> 0
    if global_flag:
        for row in range(1, len(seq_x)+1):
            alignment_matrix[row][0] = alignment_matrix[row-1][0] + scoring_matrix[seq_x[row-1]]['-']
        
        for col in range(1, len(seq_y)+1):
            alignment_matrix[0][col] = alignment_matrix[0][col-1] + scoring_matrix['-'][seq_y[col-1]]
    
    for row in range(1, len(seq_x)+1):
        for col in range(1, len(seq_y)+1):
            upper_left = alignment_matrix[row-1][col-1] + scoring_matrix[seq_x[row-1]][seq_y[col-1]]
            upper = alignment_matrix[row-1][col] + scoring_matrix[seq_x[row-1]]['-']
            left = alignment_matrix[row][col-1] + scoring_matrix['-'][seq_y[col-1]]
            if not global_flag:
                if upper_left < 0:
                    upper_left = 0
                if upper < 0:
                    upper = 0
                if left < 0:
                    left = 0
            alignment_matrix[row][col] = max(upper_left, upper, left)
        
    return alignment_matrix

#print compute_alignment_matrix('ac', 'tag', s_matrix, True)
#print compute_alignment_matrix('aa', 'taat', s_matrix, False)
#print compute_alignment_matrix('ACTACT', 'GGACTGCTTCTGG', {'A': {'A': 2, 'C': 1, '-': 0, 'T': 1, 'G': 1}, 'C': {'A': 1, 'C': 2, '-': 0, 'T': 1, 'G': 1}, '-': {'A': 0, 'C': 0, '-': 0, 'T': 0, 'G': 0}, 'T': {'A': 1, 'C': 1, '-': 0, 'T': 2, 'G': 1}, 'G': {'A': 1, 'C': 1, '-': 0, 'T': 1, 'G': 2}}, True)
#print compute_alignment_matrix('MQNSHSGVNQLGGVFVNGRPLPDSTRQKIVELAHSGARPCDISRILQVSNGCVSKILGRYYETGSIRPRAIGGSKPRVATPEVVSKIAQYKRECPSIFAWEIRDRLLSEGVCTNDNIPSVSSINRVLRNLASEKQQMGADGMYDKLRMLNGQTGSWGTRPGWYPGTSVPGQPTQDGCQQQEGGGENTNSISSNGEDSDEAQMRLQLKRKLQRNRTSFTQEQIEALEKEFERTHYPDVFARERLAAKIDLPEARIQVWFSNRRAKWRREEKLRNQRRQASNTPSHIPISSSFSTSVYQPIPQPTTPVSSFTSGSMLGRTDTALTNTYSALPPMPSFTMANNLPMQPPVPSQTSSYSCMLPTSPSVNGRSYDTYTPPHMQTHMNSQPMGTSGTTSTGLISPGVSVPVQVPGSEPDMSQYWPRLQ', 'MRNLPCLGTAGGSGLGGIAGKPSPTMEAVEASTASHPHSTSSYFATTYYHLTDDECHSGVNQLGGVFVGGRPLPDSTRQKIVELAHSGARPCDISRILQVSNGCVSKILGRYYETGSIRPRAIGGSKPRVATAEVVSKISQYKRECPSIFAWEIRDRLLQENVCTNDNIPSVSSINRVLRNLAAQKEQQSTGSGSSSTSAGNSISAKVSVSIGGNVSNVASGSRGTLSSSTDLMQTATPLNSSESGGASNSGEGSEQEAIYEKLRLLNTQHAAGPGPLEPARAAPLVGQSPNHLGTRSSHPQLVHGNHQALQQHQQQSWPPRHYSGSWYPTSLSEIPISSAPNIASVTAYASGPSLAHSLSPPNDIESLASIGHQRNCPVATEDIHLKKELDGHQSDETGSGEGENSNGGASNIGNTEDDQARLILKRKLQRNRTSFTNDQIDSLEKEFERTHYPDVFARERLAGKIGLPEARIQVWFSNRRAKWRREEKLRNQRRTPNSTGASATSSSTSATASLTDSPNSLSACSSLLSGSAGGPSVSTINGLSSPSTLSTNVNAPTLGAGIDSSESPTPIPHIRPSCTSDNDNGRQSEDCRRVCSPCPLGVGGHQNTHHIQSNGHAQGHALVPAISPRLNFNSGSFGAMYSNMHHTALSMSDSYGAVTPIPSFNHSAVGPLAPPSPIPQQGDLTPSSLYPCHMTLRPPPMAPAHHHIVPGDGGRPAGVGLGSGQSANLGASCSGSGYEVLSAYALPPPPMASSSAADSSFSAASSASANVTPHHTIAQESCPSPCSSASHFGVAHSSGFSSDPISPAVSSYAHMSYNYASSANTMTPSSASGTSAHVAPGKQQFFASCFYSPWV', {'-': {'-': -100, 'A': -5, 'C': -5, 'B': -5, 'E': -5, 'D': -5, 'G': -5, 'F': -5, 'I': -5, 'H': -5, 'K': -5, 'M': -5, 'L': -5, 'N': -5, 'Q': -5, 'P': -5, 'S': -5, 'R': -5, 'T': -5, 'W': -5, 'V': -5, 'Y': -5, 'X': -5, 'Z': -5}, 'A': {'-': -5, 'A': 5, 'C': -5, 'B': -2, 'E': -1, 'D': -2, 'G': -1, 'F': -7, 'I': -3, 'H': -5, 'K': -5, 'M': -4, 'L': -5, 'N': -2, 'Q': -3, 'P': 0, 'S': 0, 'R': -5, 'T': 0, 'W': -11, 'V': -1, 'Y': -6, 'X': -2, 'Z': -2}, 'C': {'-': -5, 'A': -5, 'C': 9, 'B': -9, 'E': -11, 'D': -11, 'G': -7, 'F': -10, 'I': -5, 'H': -6, 'K': -11, 'M': -11, 'L': -12, 'N': -8, 'Q': -11, 'P': -6, 'S': -2, 'R': -6, 'T': -6, 'W': -13, 'V': -5, 'Y': -3, 'X': -7, 'Z': -11}, 'B': {'-': -5, 'A': -2, 'C': -9, 'B': 5, 'E': 2, 'D': 6, 'G': -2, 'F': -9, 'I': -5, 'H': 0, 'K': -1, 'M': -7, 'L': -7, 'N': 5, 'Q': -2, 'P': -5, 'S': -1, 'R': -5, 'T': -2, 'W': -8, 'V': -6, 'Y': -5, 'X': -3, 'Z': 1}, 'E': {'-': -5, 'A': -1, 'C': -11, 'B': 2, 'E': 7, 'D': 3, 'G': -3, 'F': -11, 'I': -4, 'H': -3, 'K': -3, 'M': -5, 'L': -7, 'N': -1, 'Q': 2, 'P': -4, 'S': -3, 'R': -7, 'T': -4, 'W': -13, 'V': -5, 'Y': -7, 'X': -3, 'Z': 6}, 'D': {'-': -5, 'A': -2, 'C': -11, 'B': 6, 'E': 3, 'D': 7, 'G': -2, 'F': -12, 'I': -6, 'H': -2, 'K': -3, 'M': -8, 'L': -10, 'N': 2, 'Q': -1, 'P': -6, 'S': -2, 'R': -7, 'T': -3, 'W': -12, 'V': -6, 'Y': -9, 'X': -4, 'Z': 2}, 'G': {'-': -5, 'A': -1, 'C': -7, 'B': -2, 'E': -3, 'D': -2, 'G': 6, 'F': -8, 'I': -8, 'H': -7, 'K': -6, 'M': -7, 'L': -9, 'N': -2, 'Q': -5, 'P': -4, 'S': -1, 'R': -7, 'T': -4, 'W': -12, 'V': -4, 'Y': -11, 'X': -4, 'Z': -4}, 'F': {'-': -5, 'A': -7, 'C': -10, 'B': -9, 'E': -11, 'D': -12, 'G': -8, 'F': 9, 'I': -1, 'H': -5, 'K': -11, 'M': -3, 'L': -1, 'N': -7, 'Q': -10, 'P': -8, 'S': -5, 'R': -8, 'T': -7, 'W': -3, 'V': -6, 'Y': 3, 'X': -6, 'Z': -11}, 'I': {'-': -5, 'A': -3, 'C': -5, 'B': -5, 'E': -4, 'D': -6, 'G': -8, 'F': -1, 'I': 8, 'H': -7, 'K': -5, 'M': 0, 'L': 0, 'N': -4, 'Q': -6, 'P': -7, 'S': -5, 'R': -4, 'T': -1, 'W': -11, 'V': 3, 'Y': -5, 'X': -3, 'Z': -5}, 'H': {'-': -5, 'A': -5, 'C': -6, 'B': 0, 'E': -3, 'D': -2, 'G': -7, 'F': -5, 'I': -7, 'H': 9, 'K': -4, 'M': -8, 'L': -5, 'N': 1, 'Q': 2, 'P': -3, 'S': -4, 'R': 0, 'T': -5, 'W': -6, 'V': -5, 'Y': -2, 'X': -4, 'Z': 0}, 'K': {'-': -5, 'A': -5, 'C': -11, 'B': -1, 'E': -3, 'D': -3, 'G': -6, 'F': -11, 'I': -5, 'H': -4, 'K': 6, 'M': -1, 'L': -6, 'N': 0, 'Q': -2, 'P': -5, 'S': -3, 'R': 1, 'T': -2, 'W': -9, 'V': -7, 'Y': -8, 'X': -4, 'Z': -2}, 'M': {'-': -5, 'A': -4, 'C': -11, 'B': -7, 'E': -5, 'D': -8, 'G': -7, 'F': -3, 'I': 0, 'H': -8, 'K': -1, 'M': 10, 'L': 2, 'N': -6, 'Q': -3, 'P': -6, 'S': -4, 'R': -3, 'T': -3, 'W': -10, 'V': 0, 'Y': -8, 'X': -4, 'Z': -4}, 'L': {'-': -5, 'A': -5, 'C': -12, 'B': -7, 'E': -7, 'D': -10, 'G': -9, 'F': -1, 'I': 0, 'H': -5, 'K': -6, 'M': 2, 'L': 6, 'N': -6, 'Q': -4, 'P': -6, 'S': -7, 'R': -7, 'T': -5, 'W': -5, 'V': -1, 'Y': -5, 'X': -5, 'Z': -5}, 'N': {'-': -5, 'A': -2, 'C': -8, 'B': 5, 'E': -1, 'D': 2, 'G': -2, 'F': -7, 'I': -4, 'H': 1, 'K': 0, 'M': -6, 'L': -6, 'N': 7, 'Q': -2, 'P': -4, 'S': 1, 'R': -4, 'T': -1, 'W': -7, 'V': -6, 'Y': -3, 'X': -2, 'Z': -1}, 'Q': {'-': -5, 'A': -3, 'C': -11, 'B': -2, 'E': 2, 'D': -1, 'G': -5, 'F': -10, 'I': -6, 'H': 2, 'K': -2, 'M': -3, 'L': -4, 'N': -2, 'Q': 8, 'P': -2, 'S': -4, 'R': 0, 'T': -4, 'W': -10, 'V': -5, 'Y': -9, 'X': -3, 'Z': 6}, 'P': {'-': -5, 'A': 0, 'C': -6, 'B': -5, 'E': -4, 'D': -6, 'G': -4, 'F': -8, 'I': -7, 'H': -3, 'K': -5, 'M': -6, 'L': -6, 'N': -4, 'Q': -2, 'P': 8, 'S': -1, 'R': -3, 'T': -3, 'W': -11, 'V': -4, 'Y': -11, 'X': -4, 'Z': -3}, 'S': {'-': -5, 'A': 0, 'C': -2, 'B': -1, 'E': -3, 'D': -2, 'G': -1, 'F': -5, 'I': -5, 'H': -4, 'K': -3, 'M': -4, 'L': -7, 'N': 1, 'Q': -4, 'P': -1, 'S': 6, 'R': -2, 'T': 1, 'W': -4, 'V': -4, 'Y': -5, 'X': -2, 'Z': -3}, 'R': {'-': -5, 'A': -5, 'C': -6, 'B': -5, 'E': -7, 'D': -7, 'G': -7, 'F': -8, 'I': -4, 'H': 0, 'K': 1, 'M': -3, 'L': -7, 'N': -4, 'Q': 0, 'P': -3, 'S': -2, 'R': 8, 'T': -5, 'W': -1, 'V': -6, 'Y': -8, 'X': -4, 'Z': -2}, 'T': {'-': -5, 'A': 0, 'C': -6, 'B': -2, 'E': -4, 'D': -3, 'G': -4, 'F': -7, 'I': -1, 'H': -5, 'K': -2, 'M': -3, 'L': -5, 'N': -1, 'Q': -4, 'P': -3, 'S': 1, 'R': -5, 'T': 6, 'W': -10, 'V': -2, 'Y': -5, 'X': -2, 'Z': -4}, 'W': {'-': -5, 'A': -11, 'C': -13, 'B': -8, 'E': -13, 'D': -12, 'G': -12, 'F': -3, 'I': -11, 'H': -6, 'K': -9, 'M': -10, 'L': -5, 'N': -7, 'Q': -10, 'P': -11, 'S': -4, 'R': -1, 'T': -10, 'W': 13, 'V': -12, 'Y': -4, 'X': -9, 'Z': -11}, 'V': {'-': -5, 'A': -1, 'C': -5, 'B': -6, 'E': -5, 'D': -6, 'G': -4, 'F': -6, 'I': 3, 'H': -5, 'K': -7, 'M': 0, 'L': -1, 'N': -6, 'Q': -5, 'P': -4, 'S': -4, 'R': -6, 'T': -2, 'W': -12, 'V': 7, 'Y': -6, 'X': -3, 'Z': -5}, 'Y': {'-': -5, 'A': -6, 'C': -3, 'B': -5, 'E': -7, 'D': -9, 'G': -11, 'F': 3, 'I': -5, 'H': -2, 'K': -8, 'M': -8, 'L': -5, 'N': -3, 'Q': -9, 'P': -11, 'S': -5, 'R': -8, 'T': -5, 'W': -4, 'V': -6, 'Y': 9, 'X': -6, 'Z': -8}, 'X': {'-': -5, 'A': -2, 'C': -7, 'B': -3, 'E': -3, 'D': -4, 'G': -4, 'F': -6, 'I': -3, 'H': -4, 'K': -4, 'M': -4, 'L': -5, 'N': -2, 'Q': -3, 'P': -4, 'S': -2, 'R': -4, 'T': -2, 'W': -9, 'V': -3, 'Y': -6, 'X': -4, 'Z': -3}, 'Z': {'-': -5, 'A': -2, 'C': -11, 'B': 1, 'E': 6, 'D': 2, 'G': -4, 'F': -11, 'I': -5, 'H': 0, 'K': -2, 'M': -4, 'L': -5, 'N': -1, 'Q': 6, 'P': -3, 'S': -3, 'R': -2, 'T': -4, 'W': -11, 'V': -5, 'Y': -8, 'X': -3, 'Z': 6}}, True)

def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Input:  two strings seq_x and seq_y, a dictionary of dictionaries \
    scoring_matrix, and a list of lists alignment_matrix.
    Output:  a tuple of the form (score, align_x, align_y) where score \
    is the score of the global alignment align_x and align_y.
    """
    #initialize i and j to be len(x) and len(y)
    #initialize empty strings for X and Y alignments
    seq_x_index = len(seq_x)
    seq_y_index = len(seq_y)
    align_x = ""
    align_y = ""
    
    #remember that the alignment matrix has len(x) + 1 and len(y) + 1 rows and cols
    #so we can start at S[i,j] and not S[i-1,j-1]
    while seq_x_index != 0 and seq_y_index != 0:
        if alignment_matrix[seq_x_index][seq_y_index] == alignment_matrix[seq_x_index-1][seq_y_index-1] + scoring_matrix[seq_x[seq_x_index-1]][seq_y[seq_y_index-1]]:
            align_x = seq_x[seq_x_index-1] + align_x
            align_y = seq_y[seq_y_index-1] + align_y
            seq_x_index -= 1
            seq_y_index -= 1
        else:
            if alignment_matrix[seq_x_index][seq_y_index] == alignment_matrix[seq_x_index-1][seq_y_index] + scoring_matrix[seq_x[seq_x_index-1]]['-']:
                align_x = seq_x[seq_x_index-1] + align_x
                align_y = '-' + align_y
                seq_x_index -= 1
            else:
                align_x = '-' + align_x
                align_y = seq_y[seq_y_index-1] + align_y
                seq_y_index -= 1
    
    while seq_x_index != 0:
        align_x = seq_x[seq_x_index-1] + align_x
        align_y = '-' + align_y
        seq_x_index -= 1
    
    while seq_y_index != 0:
        align_x = '-' + align_x
        align_y = seq_y[seq_y_index-1] + align_y
        seq_y_index -= 1
    
    return (alignment_matrix[len(seq_x)][len(seq_y)], align_x, align_y)

#scoring_matrix = build_scoring_matrix(set(['a', 't', 'c', 'g']), 10, 4, -6)
#alignment_matrix = compute_alignment_matrix('ac', 'tag', scoring_matrix, True)
#print compute_global_alignment('ac', 'tag', scoring_matrix, alignment_matrix)

def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Input:  two sequences seq_x and seq_y, a dictionary of dictionaries \
    scoring_matrix, and a list of lists alignment matrix.
    Output:  a tuple of the form (score, align_x, align_y) where score \
    is the score of the optimal local alignment align_x and align_y.
    """
    #this should be mostly the same as compute_global_alignment except that
    #1 - it has to search alignment_matrix for the max_val, coords of max_val first
    #I think I'll implement a helper function for this
    #2 - it has to stop computing when the score is 0
    max_val, seq_x_index, seq_y_index = compute_max_val(alignment_matrix)
    
    align_x = ""
    align_y = ""
    
    while True:
        score = alignment_matrix[seq_x_index][seq_y_index]
        if score == 0:
            break
        if score == alignment_matrix[seq_x_index-1][seq_y_index-1] + scoring_matrix[seq_x[seq_x_index-1]][seq_y[seq_y_index-1]]:
            align_x = seq_x[seq_x_index-1] + align_x
            align_y = seq_y[seq_y_index-1] + align_y
            seq_x_index -= 1
            seq_y_index -= 1
        else:
            if score == alignment_matrix[seq_x_index-1][seq_y_index] + scoring_matrix[seq_x[seq_x_index-1]]['-']:
                align_x = seq_x[seq_x_index-1] + align_x
                align_y = '-' + align_y
                seq_x_index -= 1
            else:
                align_x = '-' + align_x
                align_y = seq_y[seq_y_index-1] + align_y
                seq_y_index -= 1
    
    return (max_val, align_x, align_y)

def compute_max_val(alignment_matrix):
    """
    Helper function for compute_local_alignment. Takes an alignment matrix \
    (a list of lists), performs binary search on each list (row), and returns \
    the maximum value found and the row, col location of that value as a tuple: \
    (max_val, row, col).
    """
    max_val = float('-inf')
    max_row, max_col = -1, -1
    for row in alignment_matrix:
        best_in_row = max(row)
        col = row.index(best_in_row)
        if best_in_row > max_val:
            max_val = best_in_row
            max_row = alignment_matrix.index(row)
            max_col = col
    return (max_val, max_row, max_col)

#s_matrix = build_scoring_matrix(set(['a', 't']), 10, 4, -6)
#a_matrix = compute_alignment_matrix('aa', 'taat', s_matrix, False)
#print compute_local_alignment('aa', 'taat', s_matrix, a_matrix)