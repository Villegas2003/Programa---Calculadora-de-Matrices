from fractions import Fraction

class Matriz:
    """Matriz clase y manejo de operaciones"""

    def __init__(self, filas, columnas, datos):
        self.filas = filas
        self.columnas = columnas
        self.datos = [[Fraction(x) for x in fila] for fila in datos]

    def __str__(self):
        """Convierte matriz a string para mostrarla"""
        return '\n'.join([' '.join(map(str, fila)) for fila in self.datos])

    def eliminacion_gauss_jordan(self):
        """Realiza la eliminación de Gauss-Jordan sobre la matriz."""
        n = self.filas
        m = self.columnas
        matriz = [fila[:] for fila in self.datos]  # Copia de la matriz
        pasos = []

        for i in range(n):
            if matriz[i][i] == 0:
                for j in range(i + 1, n):
                    if matriz[j][i] != 0:
                        matriz[i], matriz[j] = matriz[j], matriz[i]
                        pasos.append(f"Intercambiar fila {i + 1} con fila {j + 1}")
                        break
                else:
                    raise ValueError("No se puede aplicar Gauss-Jordan, pivote cero sin intercambios posibles.")

            pivote = matriz[i][i]
            matriz[i] = [x / pivote for x in matriz[i]]
            pasos.append(f"Dividir fila {i + 1} por {pivote}")

            for j in range(n):
                if i != j:
                    factor = matriz[j][i]
                    matriz[j] = [matriz[j][k] - factor * matriz[i][k] for k in range(m)]
                    pasos.append(f"Restar {factor} veces la fila {i + 1} de la fila {j + 1}")

        pasos.append("Resultado final:")
        pasos.append('\n'.join([' '.join(map(str, fila)) for fila in matriz]))
        print("\n".join(pasos))
        return Matriz(n, m, matriz)

    def obtener_determinante(self):
        """Obtiene el determinante de la matriz si es cuadrada y de tamaño válido"""
        if self.filas != self.columnas:
            raise ValueError("El determinante solo se puede calcular para matrices cuadradas")

        if self.filas == 1:
            return self.datos[0][0]

        try:
            def submatriz(matriz, fila, columna):
                return [[matriz[i][j] for j in range(len(matriz)) if j != columna] for i in range(len(matriz)) if i != fila]

            return sum((-1) ** c * self.datos[0][c] * Matriz(self.filas - 1, self.columnas - 1, submatriz(self.datos, 0, c)).obtener_determinante() for c in range(self.columnas))

        except Exception as e:
            raise ValueError(f"Error al calcular el determinante: {e}")

    def obtener_inversa(self, metodo="adjuncion"):
        """Obtiene la inversa de la matriz si es cuadrada y su determinante no es cero"""
        if self.filas != self.columnas:
            raise ValueError("La inversa solo se puede calcular para matrices cuadradas")

        determinante = self.obtener_determinante()
        if determinante == 0:
            raise ValueError("La matriz no tiene inversa porque su determinante es 0")

        pasos = []

        if metodo == "adjuncion":
            try:
                n = self.filas
                adjunta = [[(-1) ** (i + j) * Matriz(n - 1, n - 1, [fila[:j] + fila[j + 1:] for fila in
                                                                    (self.datos[:i] + self.datos[i + 1:])]).obtener_determinante()
                            for j in range(n)] for i in range(n)]
                print("Paso 1: Calcular la matriz adjunta.")
                print('\n'.join([' '.join(map(str, fila)) for fila in adjunta]))
                print("Paso 2: Dividir cada elemento de la adjunta por el determinante.")
                inversa = [[Fraction(adjunta[j][i], determinante) for j in range(n)] for i in range(n)]
                return Matriz(n, n, inversa)
            except Exception as e:
                raise ValueError(f"Error al calcular la inversa por adjunción: {e}")

        elif metodo == "gauss-jordan":
            try:
                n = self.filas
                matriz_extendida = [fila + [Fraction(1) if i == j else Fraction(0) for j in range(n)] for i, fila in enumerate(self.datos)]

                pasos.append("Paso 1: Formar la matriz aumentada con la matriz identidad.")
                pasos.append('\n'.join([' '.join(map(str, fila[:n])) + " | " + ' '.join(map(str, fila[n:])) for fila in matriz_extendida]))

                for i in range(n):
                    if matriz_extendida[i][i] == 0:
                        raise ValueError(
                            "No se puede calcular la inversa con el método de Gauss-Jordan debido a un pivote cero")

                    pivote = matriz_extendida[i][i]
                    matriz_extendida[i] = [round(x / pivote, 3) for x in matriz_extendida[i]]
                    pasos.append(f"Paso {i + 2}: Hacer el pivote {pivote} igual a 1 dividiendo la fila {i + 1}.")
                    pasos.append('\n'.join([' '.join(map(str, fila[:n])) + " | " + ' '.join(map(str, fila[n:])) for fila in matriz_extendida]))

                    for j in range(n):
                        if i != j:
                            factor = matriz_extendida[j][i]
                            matriz_extendida[j] = [round(matriz_extendida[j][k] - factor * matriz_extendida[i][k], 3)
                                                   for k in range(2 * n)]
                            pasos.append(f"Reducir la fila {j + 1} usando la fila {i + 1}.")
                            pasos.append('\n'.join([' '.join(map(str, fila[:n])) + " | " + ' '.join(map(str, fila[n:])) for fila in matriz_extendida]))

                inversa = [fila[n:] for fila in matriz_extendida]
                pasos.append("Resultado final:")
                pasos.append('\n'.join([' '.join(map(str, fila)) for fila in inversa]))
                print("\n".join(pasos))
                return Matriz(n, n, inversa)
            except Exception as e:
                raise ValueError(f"Error al calcular la inversa por Gauss-Jordan: {e}")

        else:
            raise ValueError("Método no reconocido. Use 'adjuncion' o 'gauss-jordan'")

    def suma(self, other):
        """Suma de matrices"""
        if self.filas != other.filas or self.columnas != other.columnas:
            raise ValueError("Matrices deben tener las mismas dimensiones")
        resultado = [
            [self.datos[i][j] + other.datos[i][j] for j in range(self.columnas)]
            for i in range(self.filas)
        ]
        return Matriz(self.filas, self.columnas, resultado)

    def resta(self, other):
        """Resta de matrices"""
        if self.filas != other.filas or self.columnas != other.columnas:
            raise ValueError("Matrices deben tener las mismas dimensiones")
        resultado = [
            [self.datos[i][j] - other.datos[i][j] for j in range(self.columnas)]
            for i in range(self.filas)
        ]
        return Matriz(self.filas, self.columnas, resultado)

    def multiplicar_matriz(self, other):
        """Multiplicación de Matrices"""
        if self.columnas != other.filas:
            raise ValueError("El # de columnas de la primera matriz deben == al numero de filas de la segunda matriz")
        resultado = []
        explicacion = []
        for i in range(self.filas):
            fila = []
            for j in range(other.columnas):
                total = 0
                pasos = []
                for k in range(self.columnas):
                    producto = self.datos[i][k] * other.datos[k][j]
                    total += producto
                    pasos.append(f"({self.datos[i][k]}×{other.datos[k][j]})")
                fila.append(total)
                explicacion.append(f"Elemento [{i + 1},{j + 1}]: {' + '.join(pasos)} = {total}")
            resultado.append(fila)
        return Matriz(self.filas, other.columnas, resultado), explicacion

    def multiplicacion_escalar(self, escalar):
        """Multiplicación Escalar"""
        resultado = [
            [element * escalar for element in fila]
            for fila in self.datos
        ]
        return Matriz(self.filas, self.columnas, resultado)

def input_matriz():
    """Maneja el ingreso de la matriz"""
    while True:
        try:
            filas = int(input("Ingrese el # de filas: "))
            columnas = int(input("Ingrese el # de columnas: "))
            if filas <= 0 or columnas <= 0:
                print("Las dimensiones deben ser un # entero positivo")
                continue

            matriz = []
            for i in range(filas):
                while True:
                    fila = input(f"Ingrese los #'s de la fila {i + 1} (separados por un espacio): ").split()
                    if len(fila) != columnas:
                        print(f"Se esperaba {columnas} de valores, se obtuvo {len(fila)}")
                        continue
                    try:
                        matriz.append([float(x) for x in fila])
                        break
                    except ValueError:
                        print("Formato invalido")
            return Matriz(filas, columnas, matriz)
        except ValueError:
            print("Porfavor ingrese #'s validos para las dimensiones")

def menu_matrices():
    matrices = []  # Lista de matrices almacenadas
    while True:
        print("\nMenú de Matrices")
        print("1. Crear nueva matriz")
        print("2. Ver matrices actuales")
        print("3. Realizar operación")
        print("4. Salir")
        
        opcion = input("Seleccione una opción: ")
        
        if opcion == '4':
            break

        if opcion == '1':
            nueva_matriz = input_matriz()
            matrices.append(nueva_matriz)
        
        elif opcion == '2':
            mostrar_matrices(matrices)
        
        elif opcion == '3':
            realizar_operacion(matrices)
        
        else:
            print("Opción inválida, intente de nuevo.")

def mostrar_matrices(matrices):
    """Muestra las matrices actuales almacenadas"""
    if matrices:
        print("Matrices actuales:")
        for idx, matriz in enumerate(matrices, 1):
            print(f"Matriz {idx}:\n{matriz}")
    else:
        print("No hay matrices actuales.")

def realizar_operacion(matrices):
    """Realiza operaciones según la opción seleccionada"""
    if len(matrices) < 2:
        print("Necesitas al menos dos matrices para realizar una operación.")
        return

    print("\nOperaciones disponibles:")
    print("a. Sumar dos matrices")
    print("b. Restar dos matrices")
    print("c. Multiplicar dos matrices")
    print("d. Multiplicación escalar")
    print("e. Determinante")
    print("f. Inversa (Adjunción)")
    print("g. Inversa (Gauss-Jordan)")
    print("h. Eliminación Gauss-Jordan")
    print("i. Regla de Cramer")

    operacion = input("Seleccione una operación: ").strip().lower()

    if operacion == 'a':
        matriz1 = matrices[int(input("Seleccione la primera matriz (índice): ")) - 1]
        matriz2 = matrices[int(input("Seleccione la segunda matriz (índice): ")) - 1]
        resultado = matriz1.suma(matriz2)
        print("Resultado de la suma:", resultado)

    elif operacion == 'b':
        matriz1 = matrices[int(input("Seleccione la primera matriz (índice): ")) - 1]
        matriz2 = matrices[int(input("Seleccione la segunda matriz (índice): ")) - 1]
        resultado = matriz1.resta(matriz2)
        print("Resultado de la resta:", resultado)

    elif operacion == 'c':
        matriz1 = matrices[int(input("Seleccione la primera matriz (índice): ")) - 1]
        matriz2 = matrices[int(input("Seleccione la segunda matriz (índice): ")) - 1]
        resultado, explicacion = matriz1.multiplicar_matriz(matriz2)
        print("Resultado de la multiplicación:", resultado)
        print("\nExplicación paso a paso:")
        for paso in explicacion:
            print(paso)

    elif operacion == 'd':
        escalar = get_escalar()
        matriz = matrices[int(input("Seleccione la matriz (índice): ")) - 1]
        resultado = matriz.multiplicacion_escalar(escalar)
        print("Resultado de la multiplicación escalar:", resultado)

    elif operacion == 'e':
        matriz = matrices[int(input("Seleccione la matriz (índice): ")) - 1]
        determinante = matriz.obtener_determinante()
        print(f"El determinante de la matriz es: {determinante}")

    elif operacion == 'f':
        matriz = matrices[int(input("Seleccione la matriz (índice): ")) - 1]
        inversa = matriz.obtener_inversa("adjuncion")
        print("La inversa por adjunción es:", inversa)

    elif operacion == 'g':
        matriz = matrices[int(input("Seleccione la matriz (índice): ")) - 1]
        inversa = matriz.obtener_inversa("gauss-jordan")
        print("La inversa por Gauss-Jordan es:", inversa)

    elif operacion == 'h':
        matriz = matrices[int(input("Seleccione la matriz (índice): ")) - 1]
        resultado = matriz.eliminacion_gauss_jordan()
        print("Resultado de la eliminación Gauss-Jordan:", resultado)

    elif operacion == 'i':
        # Aquí agregas la implementación de la Regla de Cramer
        pass

    else:
        print("Operación no válida.")

if __name__ == "__main__":
    menu_matrices()
