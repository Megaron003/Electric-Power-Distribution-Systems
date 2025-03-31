import cmath
import math
from matplotlib import pyplot as plt

# Função para converter um número complexo para a forma fasorial (magnitude < ângulo)
def to_phasor(complex_num):
    magnitude = abs(complex_num)  # Calcula a magnitude
    angle = math.degrees(cmath.phase(complex_num))  # Calcula o ângulo em graus
    return magnitude, angle

# Função para aceitar números complexos na forma retangular ou fasorial
def input_complex(prompt):
    value = input(prompt)
    try:
        # Tenta converter para número complexo na forma retangular
        return complex(value)
    except ValueError:
        # Se falhar, tenta converter da forma fasorial (magnitude < ângulo)
        magnitude, angle = map(float, value.split('<'))
        return cmath.rect(magnitude, math.radians(angle))

# Passo 1: insira os valores oriundos do sistema
v2 = input_complex("Insira um valor para V2 (na forma retangular a+bj ou fasorial magnitude<ângulo): ")
s2 = input_complex("Insira o valor estipulado da tensão em s2 (na forma retangular a+bj ou fasorial magnitude<ângulo): ")
z12 = input_complex("Insira o valor da impedância de linha do sistema, por favor: ")

# Passo 2: Cálculo de corrente
def current_calc(s2, v2):
    I = s2 / v2
    return I

# Passo 3: Cálculo de tensão na barra de carga
def charge_bar_calc(v2, z12, I):
    v1 = v2 + z12 * I
    return v1

# Passo 4: Exibe os resultados nas formas retangular e fasorial
def display_result(complex_num, label):
    print(f"{label} (Retangular): {complex_num}")
    magnitude, angle = to_phasor(complex_num)
    print(f"{label} (Fasorial): {magnitude:.4f} < {angle:.4f}°\n")

# Configuração do método iterativo
tolerance = 1e-4  # Tolerância para convergência (precisão de 4 casas decimais)
max_iterations = 100  # Número máximo de iterações para evitar loops infinitos
iteration = 0  # Contador de iterações
converged = False  # Flag para verificar se houve convergência

# Valores iniciais
I_prev = 0  # Valor inicial da corrente (pode ser qualquer valor)
V_prev = 0  # Valor inicial da tensão (pode ser qualquer valor)

# Listas para armazenar os valores de corrente e tensão
current_values = []
voltage_values = []

while not converged and iteration < max_iterations:
    iteration += 1
    print(f"\nIteração {iteration}:")

    # Calcula a corrente
    I = current_calc(s2, v2)

    # Calcula a tensão na barra de carga
    V = charge_bar_calc(v2, z12, I)

    # Armazena os valores de corrente e tensão nas listas
    current_values.append(I)
    voltage_values.append(V)

    # Exibe os resultados
    display_result(I, "Corrente I")
    display_result(V, "Tensão V1")

    # Verifica a convergência (diferença entre valores consecutivos)
    if abs(I - I_prev) < tolerance and abs(V - V_prev) < tolerance:
        converged = True
        print("Convergência alcançada!")
    else:
        # Atualiza os valores anteriores para a próxima iteração
        I_prev = I
        V_prev = V
        # Atualiza v2 para o novo valor de tensão calculado
        v2 = V

# Passo 5: Salva os valores de corrente e tensão em um arquivo .txt
with open("resultados.txt", "w") as file:
    file.write("Corrente (I)\tTensão (V)\n")  # Cabeçalho
    for I, V in zip(current_values, voltage_values):
        file.write(f"{I}\t{V}\n")

print("Resultados salvos em 'resultados.txt'.")

# Plotagem dos fasores da corrente e da tensão na barra de carga
plt.figure()
plt.polar([0, cmath.phase(I)], [0, abs(I)], marker='o', label='Corrente I')
plt.polar([0, cmath.phase(V)], [0, abs(V)], marker='o', label='Tensão V1')
plt.title('Fasores da Corrente e Tensão na Barra de Carga')
plt.legend()
plt.show()