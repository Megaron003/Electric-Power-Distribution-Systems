import math
import cmath
import numpy as np


def to_phasor(complex_num):
    magnitude = abs(complex_num)
    angle = math.degrees(cmath.phase(complex_num))
    return magnitude, angle


def input_complex(prompt):
    value = input(prompt).strip()
    try:
        value = value.replace(" ", "")
        return complex(value)
    except ValueError:
        try:
            magnitude, angle = map(float, value.split('<'))
            return cmath.rect(magnitude, math.radians(angle))
        except ValueError:
            raise ValueError("Formato inválido. Use a forma retangular (ex: 3+4j) ou fasorial (ex: 5<53.13).")


def calcular_matriz_admitancia(z_linha, z_shunt):
    Y_linha = 1 / z_linha
    Y_shunt = 1 / z_shunt

    Y = np.array([
        [Y_linha + Y_shunt, -Y_linha],
        [-Y_linha, Y_linha + Y_shunt]
    ], dtype=complex)

    print("\nMatriz de Admitância:")
    print(f"Y11 = {Y[0, 0]:.4f}")
    print(f"Y12 = {Y[0, 1]:.4f}")
    print(f"Y21 = {Y[1, 0]:.4f}")
    print(f"Y22 = {Y[1, 1]:.4f}")

    return Y


def metodo_newton(Y, V1, V2_ini, theta2_ini, P_esp, tol=1e-6, max_iter=20):
    print("\nIniciando método de Newton-Raphson...")

    theta2 = theta2_ini
    V2 = V2_ini
    iteracao = 0

    Y21 = Y[1, 0]
    G21 = Y21.real
    B21 = Y21.imag
    Y22 = Y[1, 1]
    G22 = Y22.real
    B22 = Y22.imag

    while iteracao < max_iter:
        try:
            # Cálculo da potência ativa
            P2 = V2 ** 2 * G22 + V1 * V2 * (G21 * math.cos(theta2) + B21 * math.sin(theta2))
            DeltaP = P_esp - P2

            if abs(DeltaP) < tol:
                print(f"\nConvergência alcançada na iteração {iteracao}!")
                print(f"Tensão na barra 2: {V2:.6f} pu")
                print(f"Ângulo na barra 2: {math.degrees(theta2):.6f}°")
                print(f"Potência calculada: {P2:.6f} pu")
                return V2, theta2

            # Cálculo CORRETO das derivadas parciais
            H22 = (G21 * math.sin(theta2) - B21 * math.cos(theta2))  # -∂P2/∂θ2

            # Verificação para evitar divisão por zero
            if abs(H22) < 1e-10:
                print("Aviso: Jacobiano singular, ajustando θ2 em 0.01 rad")
                theta2 += 0.01
                continue

            # Atualização CORRETA do ângulo usando o método de Newton
            Delta_theta2 = -DeltaP / H22  # Note o sinal negativo correto

            # Atualização com fator de amortecimento para estabilidade
            theta2 += Delta_theta2 * 0.8

            print(f"Iteração {iteracao}: θ2 = {math.degrees(theta2):.6f}°, P2 = {P2:.6f}, ΔP = {DeltaP:.6f}")
            iteracao += 1

        except Exception as e:
            print(f"Erro na iteração {iteracao}: {str(e)}")
            print("Reiniciando com θ2 += 0.01 rad")
            theta2 += 0.01
            continue

    print("\nAtenção: Método não convergiu após", max_iter, "iterações!")
    return V2, theta2


def main():
    print("=== Método de Newton para Fluxo de Carga em 2 Barras ===")

    # Valores de exemplo recomendados
    print("\nRecomendação para teste:")
    print("Z_linha = 0.05+0.25j (impedância série típica)")
    print("Z_shunt = 1000j (admitância shunt alta)")
    print("V1 = 1.0, V2_ini = 1.0, θ2_ini = 0.1°")
    print("P_esp entre 0.1 e 0.8 pu para melhor convergência\n")

    z_linha = input_complex("Insira a impedância de série da linha: ")
    z_shunt = input_complex("Insira a impedância shunt da linha: ")

    Y = calcular_matriz_admitancia(z_linha, z_shunt)

    V1 = float(input("\nTensão na barra 1 (pu): "))
    V2_ini = float(input("Tensão inicial na barra 2 (pu): "))
    theta2_ini = math.radians(float(input("Ângulo inicial na barra 2 (graus): ")))

    P_esp = float(input("\nPotência ativa especificada na barra 2 (pu): "))
    tol = float(input("Tolerância de convergência (padrão 1e-6): ") or "1e-6")

    V2_final, theta2_final = metodo_newton(Y, V1, V2_ini, theta2_ini, P_esp, tol)

    # Resultado final
    print("\n=== Resultado Final ===")
    print(f"Ângulo calculado: {math.degrees(theta2_final):.6f}°")
    print(f"Tensão calculada: {V2_final:.6f} pu")

    # Verificação
    P_final = V2_final ** 2 * Y[1, 1].real + V1 * V2_final * (
            Y[1, 0].real * math.cos(theta2_final) + Y[1, 0].imag * math.sin(theta2_final))
    print(f"Potência ativa resultante: {P_final:.6f} pu")
    print(f"Erro: {P_esp - P_final:.2e} pu")


if __name__ == "__main__":
    main()