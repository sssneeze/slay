#!/usr/bin/env python3

import csv
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def load_data(filename="C:/Users/asave/slay/symmetric_acceleration_results.csv"):
    methods = {}
    
    with open(filename, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            method = row['method']
            if method not in methods:
                methods[method] = {
                    'iterations': [],
                    'residuals': [],
                    'times': []
                }
            methods[method]['iterations'].append(int(row['iteration']))
            methods[method]['residuals'].append(float(row['residual']))
            methods[method]['times'].append(float(row['time_ms']))
    
    return methods


def plot_convergence_by_iterations(methods, filename='symmetric_convergence_iterations.png'):
    plt.figure(figsize=(10, 6))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(methods)))
    
    for idx, (method, data) in enumerate(methods.items()):
        plt.semilogy(data['iterations'], data['residuals'], 
                    marker='o', markersize=4, linewidth=2, 
                    label=method, color=colors[idx])
    
    plt.xlabel('Номер итерации', fontsize=12)
    plt.ylabel('Невязка (логарифмическая шкала)', fontsize=12)
    plt.title('Сходимость итерационных методов\n(симметричная матрица Пуассона, N=100)', 
             fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()


def plot_convergence_by_time(methods, filename='symmetric_convergence_time.png'):
    plt.figure(figsize=(10, 6))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(methods)))
    
    for idx, (method, data) in enumerate(methods.items()):
        plt.semilogy(data['times'], data['residuals'], 
                    marker='s', markersize=4, linewidth=2, 
                    label=method, color=colors[idx])
    
    plt.xlabel('Время, мс', fontsize=12)
    plt.ylabel('Невязка (логарифмическая шкала)', fontsize=12)
    plt.title('Сходимость итерационных методов по времени\n(симметричная матрица Пуассона, N=100)', 
             fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()



def main():
    csv_file = "C:/Users/asave/slay/symmetric_acceleration_results.csv"
    

    methods = load_data(csv_file)
    
    for method, data in methods.items():
        print(f"  - {method}: {len(data['iterations'])} точек")
    
    plot_convergence_by_iterations(methods)
    plot_convergence_by_time(methods)


if __name__ == '__main__':
    main()