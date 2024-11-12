import numpy as np

def read_file(input):
    with open(input, 'r') as f:
        seq1 = f.readline().strip()
        seq2 = f.readline().strip()
        gap = int(f.readline().strip() or -2)
        mismatch = int(f.readline().strip() or -1)
        match = int(f.readline().strip() or 1)
    return seq1, seq2, gap, mismatch, match

def create_matrix(seq1, seq2, gap, mismatch, match):
    len1 = len(seq1) + 1
    len2 = len(seq2) + 1
    matrix = np.zeros((len1, len2))
    return matrix, len1, len2

def calculate_score(matrix, seq1, seq2, gap, mismatch, match, len1, len2):
    for i in range(1, len1):
        matrix[i][0] = matrix[i-1][0] + gap

    for j in range(1, len2):
        matrix[0][j] = matrix[0][j-1] + gap

    for i in range(1, len1):
        for j in range(1, len2):
            match_mismatch = match if seq1[i-1] == seq2[j-1] else mismatch
            matrix[i][j] = max(matrix[i-1][j-1] + match_mismatch,
                               matrix[i-1][j] + gap,
                               matrix[i][j-1] + gap)

    return matrix

def print_matrix_with_headers(matrix, seq1, seq2, file):
    print("--------------------------------------------------------------------------------", file=file)
    print("** valores de score **", file=file)
    print("================================================================================", file=file)

    # Cabeçalho das colunas com o espaço inicial para o alinhamento

    

    # A matriz precisa ser invertida para seguir a ordem de cima para baixo conforme pedido
    for i, row in enumerate(matrix[::-1]):
        # Inverte o índice para manter o cabeçalho na ordem correta "C T A U"
        header = seq1[len(seq1) - 1 - i] if i < len(seq1) else "U"
        print(f"{header:<3}", end=" ", file=file)

        # Imprime cada valor da linha invertida
        for value in row:
            print(f"{int(value):^3}", end=" ", file=file)
        print(file=file)

    # Cabeçalho inferior "X U T C G"
    print("X", "U", *seq2, sep="   ", file=file)




    print("================================================================================", file=file)

def traceback(matrix, seq1, seq2, len1, len2, gap, mismatch, match):
    align1 = ""
    align2 = ""
    
    # Encontrar o maior valor na última coluna
    max_value = -float('inf')
    max_row = len1 - 1
    for i in range(len1):
        if matrix[i][len2 - 1] > max_value:
            max_value = matrix[i][len2 - 1]
            max_row = i

    # Começar o rastreamento a partir da célula com o maior valor na última coluna
    i, j = max_row, len2 - 1
    #i, j = len1 - 1, len2 - 1   #Aqui é pra caso eu queira começar no canto superior direito.

    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i-1][j] + gap:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1

    return align1, align2

def print_alignment_to_file(align1, align2, score, match, mismatch, gap, output_file):
    print(f"Score = {int(score)}", file=output_file)
    print(f"** Match = {match} | mismatch = {mismatch} | Gap = {gap}**", file=output_file)
    print("------------------------------------------------------------------", file=output_file)
    print("Alinhamento\n" , file=output_file)
    print(align1, file=output_file)
    print(align2, file=output_file)

def smith_waterman_global(input, output):
    seq1, seq2, gap, mismatch, match = read_file(input)
    matrix, len1, len2 = create_matrix(seq1, seq2, gap, mismatch, match)
    matrix = calculate_score(matrix, seq1, seq2, gap, mismatch, match, len1, len2)

    with open(output, 'w') as output_file:
        print_matrix_with_headers(matrix, seq1, seq2, output_file)

        align1, align2 = traceback(matrix, seq1, seq2, len1, len2, gap, mismatch, match)
        score = matrix[len1-1][len2-1]   #score local
        #score = np.max(matrix)   #score global 

        print_alignment_to_file(align1, align2, score, match, mismatch, gap, output_file)

input = 'input.txt'
output = 'output.txt'
smith_waterman_global(input, output)
