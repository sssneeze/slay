#!/usr/bin/env python3
import csv
import matplotlib.pyplot as plt


with open("C:/Users/asave/slay/methods_results.csv", 'r') as f:
    reader = csv.DictReader(f)
    data = list(reader)

methods = [r['method'] for r in data]
iterations = [int(r['iterations']) for r in data]
times = [float(r['time_ms']) for r in data]
residuals = [float(r['residual']) for r in data]


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.bar(methods, iterations, edgecolor='black')
ax1.set_xlabel('Метод')
ax1.set_ylabel('Итерации')
ax1.tick_params(axis='x', rotation=15)

ax2.bar(methods, times, edgecolor='black')
ax2.set_xlabel('Метод')
ax2.set_ylabel('Время, мс')
ax2.tick_params(axis='x', rotation=15)

plt.tight_layout()
plt.savefig('plot_iter_time.png', dpi=300)




plt.figure(figsize=(6, 4))
plt.bar(methods, residuals, edgecolor='black')
plt.yscale('log')
plt.xlabel('Метод')
plt.ylabel('Невязка')
plt.tick_params(axis='x', rotation=15)
plt.tight_layout()
plt.savefig('plot_residual.png', dpi=300)