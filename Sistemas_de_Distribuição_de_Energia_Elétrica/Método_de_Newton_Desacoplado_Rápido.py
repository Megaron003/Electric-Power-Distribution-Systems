import math
import cmath
import numpy as np
from scipy.misc import derivative


# Função para converter para forma fasorial
def to_phasor(complex_num):
    return abs(complex_num), math.degrees(cmath.phase(complex_num))


# Função para entrada de números complexos
def input_complex(prompt):
    value = input(prompt).strip()
    try:
        return complex(value.replace(" ", ""))
    except ValueError:
        try:
            mag, ang = map(float, value.split('<'))
            return cmath.rect(mag, math.radians(ang))
        except:
            raise ValueError("Formato inválido. Use 3+4j ou 5<53.13")


# Cálculo da matriz de admitância
def matrix_calc(z, shunt):
    Y = 1 / z
    Y_shunt = 1 / shunt
    matrix = [
        [Y + Y_shunt, -Y],
        [-Y, Y + Y_shunt]
    ]
    print("\nMatriz de Admitância:")
    for row in matrix:
        print([f"{elem.real:.4f}{elem.imag:+.4f}j" for elem in row])
    return matrix


# Função para calcular a potência ativa
def power_equation(theta, matrix, v1, v2):
    G21 = matrix[1][0].real
    B21 = matrix[1][0].imag
    G22 = matrix[1][1].real
    return (v2 ** 2) * G22 + v1 * v2 * (G21 * math.cos(theta) + B21 * math.sin(theta))


# Método de Newton com derivada numérica
def newton_method(P_spec, P_calc, epsilon, matrix, v1, v2, theta_init):
    theta = theta_init
    iteration = 0
    max_iter = 100

    while True:
        # Calcula a derivada numérica usando scipy
        dp_dtheta = derivative(
            lambda x: power_equation(x, matrix, v1, v2),
            theta,
            dx=1e-6
        )

        # Verifica se a derivada é muito pequena
        if abs(dp_dtheta) < 1e-10:
            print("Aviso: Derivada muito pequena, ajustando theta...")
            theta += 0.01
            continue

        Delta_P = P_spec - P_calc
        if abs(Delta_P) < epsilon:
            print(f"\nConvergência alcançada na iteração {iteration}!")
            print(f"Ângulo final: {math.degrees(theta):.6f}°")
            print(f"Potência calculada: {P_calc:.6f}")
            return theta

        if iteration >= max_iter:
            print("\nAtenção: Número máximo de iterações atingido!")
            return theta

        # Atualiza theta usando Newton-Raphson
        delta_theta = Delta_P / dp_dtheta
        theta += delta_theta

        # Recalcula a potência com novo theta
        P_calc = power_equation(theta, matrix, v1, v2)

        print(f"Iteração {iteration}: θ = {math.degrees(theta):.6f}°, P = {P_calc:.6f}, ΔP = {Delta_P:.6f}")
        iteration += 1


# Programa principal
if __name__ == "__main__":
    print("=== Método de Newton para Fluxo de Carga com Derivada Numérica ===")

    # Entrada de dados
    print("\nInsira a impedância de linha (ex: 0.1+0.5j ou 0.5<78.69):")
    z_line = input_complex("Z_linha = ")

    print("\nInsira a impedância shunt (ex: 1000j):")
    z_shunt = input_complex("Z_shunt = ")

    # Calcula matriz Y
    Y_matrix = matrix_calc(z_line, z_shunt)

    # Condições iniciais
    v1 = float(input("\nTensão na barra 1 (pu): "))
    v2 = float(input("Tensão na barra 2 (pu): "))
    theta = math.radians(float(input("Ângulo inicial na barra 2 (graus): ")))

    # Potência especificada
    P_spec = float(input("\nPotência ativa especificada na barra 2 (pu): "))
    tolerance = float(input("Tolerância de convergência (ex: 1e-6): ") or "1e-6")

    # Cálculo inicial
    P_initial = power_equation(theta, Y_matrix, v1, v2)
    print(f"\nPotência inicial: {P_initial:.6f}")

    # Executa o método de Newton
    final_theta = newton_method(P_spec, P_initial, tolerance, Y_matrix, v1, v2, theta)

    # Resultados finais
    final_power = power_equation(final_theta, Y_matrix, v1, v2)
    print("\n=== Resultados Finais ===")
    print(f"Ângulo calculado: {math.degrees(final_theta):.6f}°")
    print(f"Potência calculada: {final_power:.6f} pu")
    print(f"Erro: {P_spec - final_power:.2e} pu")