import math
import cmath
import numpy as np
from tkinter import *
from tkinter import ttk, messagebox


class PowerFlowGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Fluxo de Carga - Método de Newton-Raphson")
        self.root.geometry("800x600")

        # Variáveis do sistema
        self.v1 = 1.0
        self.v2 = 1.0
        self.theta2 = 0.0
        self.Y_matrix = None
        self.P_esp = -0.2
        self.Q_esp = -0.01
        self.G21 = 0.0
        self.B21 = 0.0
        self.G22 = 0.0
        self.B22 = 0.0

        self.create_widgets()

    def create_widgets(self):
        # Frame principal
        mainframe = ttk.Frame(self.root, padding="10")
        mainframe.grid(column=0, row=0, sticky=(N, W, E, S))

        # Entradas de dados
        ttk.Label(mainframe, text="Impedância de Linha (ex: 0.1+0.5j ou 0.5<78.69):").grid(column=1, row=1, sticky=W)
        self.z_line = ttk.Entry(mainframe, width=20)
        self.z_line.grid(column=2, row=1, sticky=(W, E))

        ttk.Label(mainframe, text="Impedância Shunt (ex: 1000j):").grid(column=1, row=2, sticky=W)
        self.z_shunt = ttk.Entry(mainframe, width=20)
        self.z_shunt.grid(column=2, row=2, sticky=(W, E))

        ttk.Label(mainframe, text="Tensão Inicial Barra 2 (pu):").grid(column=1, row=3, sticky=W)
        self.v2_entry = ttk.Entry(mainframe, width=20)
        self.v2_entry.grid(column=2, row=3, sticky=(W, E))
        self.v2_entry.insert(0, "1.0")

        ttk.Label(mainframe, text="Ângulo Inicial Barra 2 (rad):").grid(column=1, row=4, sticky=W)
        self.theta2_entry = ttk.Entry(mainframe, width=20)
        self.theta2_entry.grid(column=2, row=4, sticky=(W, E))
        self.theta2_entry.insert(0, "0.0")

        ttk.Label(mainframe, text="Potência Ativa Especificada (pu):").grid(column=1, row=5, sticky=W)
        self.P_esp_entry = ttk.Entry(mainframe, width=20)
        self.P_esp_entry.grid(column=2, row=5, sticky=(W, E))
        self.P_esp_entry.insert(0, "-0.2")

        ttk.Label(mainframe, text="Potência Reativa Especificada (pu):").grid(column=1, row=6, sticky=W)
        self.Q_esp_entry = ttk.Entry(mainframe, width=20)
        self.Q_esp_entry.grid(column=2, row=6, sticky=(W, E))
        self.Q_esp_entry.insert(0, "-0.01")

        # Botões
        ttk.Button(mainframe, text="Calcular", command=self.run_calculation).grid(column=2, row=7, sticky=W)
        ttk.Button(mainframe, text="Limpar", command=self.clear_fields).grid(column=2, row=7, sticky=E)

        # Área de resultados
        self.results_text = Text(mainframe, width=80, height=20, wrap=WORD)
        self.results_text.grid(column=1, row=8, columnspan=2, pady=10)

        # Barra de rolagem
        scrollbar = ttk.Scrollbar(mainframe, orient=VERTICAL, command=self.results_text.yview)
        scrollbar.grid(column=3, row=8, sticky=(N, S))
        self.results_text['yscrollcommand'] = scrollbar.set

    def input_complex(self, value):
        """Converte string para número complexo"""
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
        Y_shunt = 0 if shunt == 0 else 1 / shunt

        self.Y_matrix = np.array([
            [Y + Y_shunt, -Y],
            [-Y, Y + Y_shunt]
        ], dtype=complex)

        return self.Y_matrix

    def calculate_power(self, theta_2=None):
        """Calcula P e Q na barra 2"""
        if theta_2 is None:
            theta_2 = self.theta2

        self.G21 = self.Y_matrix[1, 0].real
        self.B21 = self.Y_matrix[1, 0].imag
        self.G22 = self.Y_matrix[1, 1].real
        self.B22 = self.Y_matrix[1, 1].imag

        P = (self.v2 ** 2) * self.G22 + self.v1 * self.v2 * (
                    self.G21 * math.cos(theta_2) + self.B21 * math.sin(theta_2))
        Q = -(self.v2 ** 2) * self.B22 + self.v1 * self.v2 * (
                    self.G21 * math.sin(theta_2) - self.B21 * math.cos(theta_2))

        return P, Q

    def calculate_jacobian(self):
        """Calcula a matriz Jacobiana"""
        self.H22 = self.v1 * self.v2 * (-self.G21 * math.sin(self.theta2) + self.B21 * math.cos(self.theta2))
        self.N22 = 2 * self.v2 * self.G22 + self.v1 * (
                    self.G21 * math.cos(self.theta2) + self.B21 * math.sin(self.theta2))
        self.M22 = self.v1 * self.v2 * (self.G21 * math.cos(self.theta2) + self.B21 * math.sin(self.theta2))
        self.L22 = -2 * self.v2 * self.B22 + self.v1 * (
                    self.G21 * math.sin(self.theta2) - self.B21 * math.cos(self.theta2))

        return np.array([[self.H22, self.N22], [self.M22, self.L22]])

    def newton_raphson(self, tol=1e-6, max_iter=20):
        """Implementa o método de Newton-Raphson"""
        self.results_text.insert(END, "\nIniciando método de Newton-Raphson...\n")

        for iteration in range(max_iter):
            P, Q = self.calculate_power()
            DeltaP = self.P_esp - P
            DeltaQ = self.Q_esp - Q

            if abs(DeltaP) < tol and abs(DeltaQ) < tol:
                self.results_text.insert(END, f"\nConvergência alcançada na iteração {iteration}!\n")
                self.results_text.insert(END, f"Tensão na barra 2: {self.v2:.6f} pu\n")
                self.results_text.insert(END, f"Ângulo na barra 2: {self.theta2:.6f} rad\n")
                return True

            J = self.calculate_jacobian()

            try:
                delta = np.linalg.solve(J, np.array([DeltaP, DeltaQ]))
            except np.linalg.LinAlgError:
                self.results_text.insert(END, "Erro: Jacobiano singular\n")
                return False

            self.theta2 += delta[0]
            self.v2 += delta[1]

            self.results_text.insert(END, f"\nIteração {iteration}:\n")
            self.results_text.insert(END, f"  Δθ = {delta[0]:.6f} rad, ΔV = {delta[1]:.6f} pu\n")
            self.results_text.insert(END, f"  P = {P:.6f} pu, Q = {Q:.6f} pu\n")
            self.results_text.insert(END, f"  ΔP = {DeltaP:.6f} pu, ΔQ = {DeltaQ:.6f} pu\n")

            self.results_text.insert(END, "\n  Elementos da Matriz Jacobiana:\n")
            self.results_text.insert(END, f"  H22 = {self.H22:.6f}\n")
            self.results_text.insert(END, f"  N22 = {self.N22:.6f}\n")
            self.results_text.insert(END, f"  M22 = {self.M22:.6f}\n")
            self.results_text.insert(END, f"  L22 = {self.L22:.6f}\n")

            self.results_text.insert(END, "\n  Coeficientes da matriz de Admitância:\n")
            self.results_text.insert(END, f"  G21 = {self.G21:.6f}\n")
            self.results_text.insert(END, f"  B21 = {self.B21:.6f}\n")
            self.results_text.insert(END, f"  G22 = {self.G22:.6f}\n")
            self.results_text.insert(END, f"  B22 = {self.B22:.6f}\n")

        self.results_text.insert(END, f"\nAtenção: Método não convergiu após {max_iter} iterações!\n")
        return False

    def run_calculation(self):
        """Executa o fluxo de carga"""
        try:
            self.results_text.delete(1.0, END)

            # Obter valores da interface
            z_line = self.input_complex(self.z_line.get())
            z_shunt = self.input_complex(self.z_shunt.get())
            self.v2 = float(self.v2_entry.get())
            self.theta2 = float(self.theta2_entry.get())
            self.P_esp = float(self.P_esp_entry.get())
            self.Q_esp = float(self.Q_esp_entry.get())

            # Calcular matriz de admitância
            Y_matrix = self.matrix_calc(z_line, z_shunt)
            self.results_text.insert(END, "\nMatriz de Admitância:\n")
            self.results_text.insert(END, str(np.array_str(Y_matrix, precision=4, suppress_small=True)) + "\n")

            # Executar método
            self.newton_raphson()

            # Resultados finais
            P_final, Q_final = self.calculate_power()
            self.results_text.insert(END, "\n=== Resultados Finais ===\n")
            self.results_text.insert(END, f"Tensão final na barra 2: {self.v2:.6f} pu\n")
            self.results_text.insert(END, f"Ângulo final na barra 2: {self.theta2:.6f} rad\n")
            self.results_text.insert(END, f"Potência ativa calculada: {P_final:.6f} pu\n")
            self.results_text.insert(END, f"Potência reativa calculada: {Q_final:.6f} pu\n")

        except Exception as e:
            messagebox.showerror("Erro", f"Ocorreu um erro:\n{str(e)}")

    def clear_fields(self):
        """Limpa todos os campos de entrada"""
        self.z_line.delete(0, END)
        self.z_shunt.delete(0, END)
        self.v2_entry.delete(0, END)
        self.v2_entry.insert(0, "1.0")
        self.theta2_entry.delete(0, END)
        self.theta2_entry.insert(0, "0.0")
        self.P_esp_entry.delete(0, END)
        self.P_esp_entry.insert(0, "-0.2")
        self.Q_esp_entry.delete(0, END)
        self.Q_esp_entry.insert(0, "-0.01")
        self.results_text.delete(1.0, END)


# Iniciar a aplicação
if __name__ == "__main__":
    root = Tk()
    app = PowerFlowGUI(root)
    root.mainloop()