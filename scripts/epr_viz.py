import matplotlib.pyplot as plt
import argparse
import numpy as np

# Настройка парсинга аргументов командной строки
parser = argparse.ArgumentParser()
parser.add_argument("file", type=str, help="Путь к файлу с данными")
args = parser.parse_args()

angles = []
values_db = []

filename = args.file 

try:
    with open(filename, "r") as file:
        for line in file:
            angle, value = map(float, line.split())
            angles.append(angle)
            values_db.append(value)
except FileNotFoundError:
    print(f"Файл '{filename}' не найден.")
    exit()

values_db = np.array(values_db)
values_db = 10*np.log10(values_db/np.pi)
# Построение графика
plt.figure(figsize=(10, 6))
plt.plot(angles, values_db, label="Magnitude (dB)", color="blue")

# Настройка графика
plt.xlabel("Угол (градусы)", fontsize=12)
plt.ylabel("дБ", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)
plt.tight_layout()

plt.show()