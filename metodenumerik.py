# Fungsi untuk menghitung determinan matriks
def determinan(matriks):
    n = len(matriks)
    if n == 1:
        return matriks[0][0]
    elif n == 2:
        return matriks[0][0] * matriks[1][1] - matriks[0][1] * matriks[1][0]
    else:
        det = 0
        for i in range(n):
            sub_matriks = [row[:i] + row[i+1:] for row in (matriks[:i] + matriks[i+1:])]
            det += ((-1) ** i) * matriks[0][i] * determinan(sub_matriks)
        return det

# Metode Matriks Balikan
def metode_balikan(matriks_A, vektor_b):
    n = len(matriks_A)
    det_A = determinan(matriks_A)
    if det_A == 0:
        print("Matriks tidak dapat dibalikan.")
        return None
    
    matriks_balikan = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            sub_matriks = [row[:j] + row[j+1:] for row in (matriks_A[:i] + matriks_A[i+1:])]
            matriks_balikan[i][j] = ((-1) ** (i + j)) * determinan(sub_matriks) / det_A
    
    solusi = [sum(matriks_balikan[i][j] * vektor_b[j] for j in range(n)) for i in range(n)]
    return solusi

# Metode Dekomposisi LU Gauss
def dekomposisi_LU_Gauss(matriks_A, vektor_b):
    n = len(matriks_A)
    matriks_L = [[0 for _ in range(n)] for _ in range(n)]
    matriks_U = [[0 for _ in range(n)] for _ in range(n)]
    
    for i in range(n):
        for j in range(i, n):
            sum = 0
            for k in range(i):
                sum += matriks_L[i][k] * matriks_U[k][j]
            matriks_U[i][j] = matriks_A[i][j] - sum
        
        for j in range(i, n):
            if matriks_U[i][i] == 0:
                print("Matriks tidak dapat difaktorisasi.")
                return None
            sum = 0
            for k in range(i):
                sum += matriks_L[j][k] * matriks_U[k][i]
            matriks_L[j][i] = (matriks_A[j][i] - sum) / matriks_U[i][i]
    
    for i in range(n):
        matriks_L[i][i] = 1
    
    vektor_y = [0 for _ in range(n)]
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += matriks_L[i][j] * vektor_y[j]
        vektor_y[i] = vektor_b[i] - sum
    
    solusi = [0 for _ in range(n)]
    for i in range(n-1, -1, -1):
        sum = 0
        for j in range(i+1, n):
            sum += matriks_U[i][j] * solusi[j]
        solusi[i] = (vektor_y[i] - sum) / matriks_U[i][i]
    
    return solusi

# Metode Dekomposisi Crout
def dekomposisi_Crout(matriks_A, vektor_b):
    n = len(matriks_A)
    matriks_L = [[0 for _ in range(n)] for _ in range(n)]
    matriks_U = [[0 for _ in range(n)] for _ in range(n)]
    
    for i in range(n):
        for j in range(i, n):
            sum = 0
            for k in range(i):
                sum += matriks_L[i][k] * matriks_U[k][j]
            matriks_U[i][j] = matriks_A[i][j] - sum
        
        for j in range(i, n):
            if matriks_U[j][j] == 0:
                print("Matriks tidak dapat difaktorisasi (elemen diagonal utama bernilai nol).")
                return None
            sum = 0
            for k in range(j):
                sum += matriks_L[i][k] * matriks_U[k][j]
            matriks_L[i][j] = (matriks_A[i][j] - sum) / matriks_U[j][j]
    
    for i in range(n):
        matriks_L[i][i] = 1
    
    vektor_y = [0 for _ in range(n)]
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += matriks_L[i][j] * vektor_y[j]
        vektor_y[i] = vektor_b[i] - sum
    
    solusi = [0 for _ in range(n)]
    for i in range(n-1, -1, -1):
        sum = 0
        for j in range(i+1, n):
            sum += matriks_U[i][j] * solusi[j]
        solusi[i] = (vektor_y[i] - sum) / matriks_U[i][i]
    
    return solusi

# Kode untuk pengujian
matriks_A = [[3, 2, 1], [1, 3, 2], [2, 1, 3]]
vektor_b = [6, 5, 7]

print("Metode Matriks Balikan:")
solusi_balikan = metode_balikan(matriks_A, vektor_b)
print(solusi_balikan)

print("\nMetode Dekomposisi LU Gauss:")
solusi_LU_Gauss = dekomposisi_LU_Gauss(matriks_A, vektor_b)
print(solusi_LU_Gauss)

print("\nMetode Dekomposisi Crout:")
solusi_Crout = dekomposisi_Crout(matriks_A, vektor_b)
print(solusi_Crout)