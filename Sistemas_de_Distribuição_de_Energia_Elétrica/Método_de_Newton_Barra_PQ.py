import math
import cmath
import numpy as np
from scipy.misc import derivative


class PowerFlowNewton:
    def __init__(self):
        self.v1 = 1.0  # Tensão na barra slack (pu)
        self.v2 = 1.0  # Tensão inicial na barra PQ (pu)
        self.theta2 = 0.0  # Ângulo inicial (rad)
        self.Y_matrix = None
        self.P_esp = -0.2
        self.Q_esp = -0.01
        # Inicializando os coeficientes
        self.G21 = 0.0
        self.B21 = 0.0
        self.G22 = 0.0
        self.B22 = 0.0

    def input_complex(self, prompt):
        """Função para entrada de números complexos"""
        value = input(prompt).strip()
        try:
            return complex(value.replace(" ", ""))
        except ValueError:
            try:
                mag, ang = map(float, value.split('<'))
                return cmath.rect(mag, math.radians(ang))
            except:
                raise ValueError("Formato inválido. Use 3+4j ou 5<53.13")

    def matrix_calc(self, z, shunt):
        """Calcula a matriz de admitância"""
        Y = 1 / z

        if shunt == 0:
            Y_shunt = 0
        else:
            Y_shunt = 1 / shunt

        self.Y_matrix = np.array([
            [Y + Y_shunt, -Y],
            [-Y, Y + Y_shunt]
        ], dtype=complex)

        print("\nMatriz de Admitância:")
        print(np.array_str(self.Y_matrix, precision=4, suppress_small=True))
        return self.Y_matrix

    def calculate_power(self, theta_2=None):
        """Calcula P e Q na barra 2"""
        if theta_2 is None:
            theta_2 = self.theta2

        self.G21 = self.Y_matrix[1, 0].real
        self.B21 = self.Y_matrix[1, 0].imag
        self.G22 = self.Y_matrix[1, 1].real
        self.B22 = self.Y_matrix[1, 1].imag

        # Potência ativa
        P = (self.v2 ** 2) * self.G22 + self.v1 * self.v2 * (self.G21 * math.cos(theta_2) + self.B21 * math.sin(theta_2))

        # Potência reativa
        Q = -(self.v2 ** 2) * self.B22 + self.v1 * self.v2 * (self.G21 * math.sin(theta_2) - self.B21 * math.cos(theta_2))

        return P, Q

    def calculate_jacobian(self):
        """Calcula a matriz Jacobiana"""
        # Elementos da Jacobiana
        self.H22 = self.v1 * self.v2 * (-self.G21 * math.sin(self.theta2) + self.B21 * math.cos(self.theta2))
        self.N22 = 2 * self.v2 * self.G22 + self.v1 * (self.G21 * math.cos(self.theta2) + self.B21 * math.sin(self.theta2))
        self.M22 = self.v1 * self.v2 * (self.G21 * math.cos(self.theta2) + self.B21 * math.sin(self.theta2))
        self.L22 = -2 * self.v2 * self.B22 + self.v1 * (self.G21 * math.sin(self.theta2) - self.B21 * math.cos(self.theta2))

        J = np.array([
            [self.H22, self.N22],
            [self.M22, self.L22]
        ])
        return J

    def newton_raphson(self, tol=1e-6, max_iter=20):
        """Implementa o método de Newton-Raphson"""
        print("\nIniciando método de Newton-Raphson...")

        for iteration in range(max_iter):
            P, Q = self.calculate_power()
            DeltaP = self.P_esp - P
            DeltaQ = self.Q_esp - Q

            # Verificar convergência
            if abs(DeltaP) < tol and abs(DeltaQ) < tol:
                print(f"\nConvergência alcançada na iteração {iteration}!")
                print(f"Tensão na barra 2: {self.v2:.6f} pu")
                print(f"Ângulo na barra 2: {self.theta2:.6f} rad")
                return True

            # Calcular Jacobiano
            J = self.calculate_jacobian()

            # Resolver sistema linear
            try:
                delta = np.linalg.solve(J, np.array([DeltaP, DeltaQ]))
            except np.linalg.LinAlgError:
                print("Erro: Jacobiano singular")
                return False

            # Atualizar variáveis
            self.theta2 += delta[0]
            self.v2 += delta[1]

            print(f"Iteração {iteration}:")
            print(f"  Δθ = {delta[0]:.6f} rad, ΔV = {delta[1]:.6f} pu")
            print(f"  P = {P:.6f} pu, Q = {Q:.6f} pu")
            print(f"  ΔP = {DeltaP:.6f} pu, ΔQ = {DeltaQ:.6f} pu")

            print("\n  Elementos da Matriz Jacobiana:")
            print(f"  H22 = {self.H22:.6f}")
            print(f"  N22 = {self.N22:.6f}")
            print(f"  M22 = {self.M22:.6f}")
            print(f"  L22 = {self.L22:.6f}")

            print("\n  Coeficientes da matriz de Admitância:")
            print(f"  G21 = {self.G21:.6f}")
            print(f"  B21 = {self.B21:.6f}")
            print(f"  G22 = {self.G22:.6f}")
            print(f"  B22 = {self.B22:.6f}\n")

        print("\nAtenção: Método não convergiu após", max_iter, "iterações!")
        return False

    def run(self):
        """Executa o fluxo de carga"""
        print("=== Método de Newton para Fluxo de Carga em Sistema PQ ===")

        # Entrada de dados
        print("\nInsira a impedância de linha (ex: 0.1+0.5j ou 0.5<78.69):")
        z_line = self.input_complex("Z_linha = ")

        print("\nInsira a impedância shunt (ex: 1000j):")
        z_shunt = self.input_complex("Z_shunt = ")

        # Calcular matriz de admitância
        self.matrix_calc(z_line, z_shunt)

        # Configurar condições iniciais
        self.v2 = float(input("\nTensão inicial na barra 2 (pu): "))
        self.theta2 = float(input("Ângulo inicial na barra 2 (radianos): "))

        # Potências especificadas
        self.P_esp = float(input("\nPotência ativa especificada na barra 2 em P.U. (P_esp): "))
        self.Q_esp = float(input("Potência reativa especificada na barra 2 em P.U. (Q_esp): "))

        # Executar método
        self.newton_raphson()

        # Resultados finais
        P_final, Q_final = self.calculate_power()
        print("\n=== Resultados Finais ===")
        print(f"Tensão final na barra 2: {self.v2:.6f} pu")
        print(f"Ângulo final na barra 2: {self.theta2:.6f} rad")
        print(f"Potência ativa calculada: {P_final:.6f} pu")
        print(f"Potência reativa calculada: {Q_final:.6f} pu")


if __name__ == "__main__":
    power_flow = PowerFlowNewton()
    power_flow.run()