import pandas as pd
import matplotlib.pyplot as plt

# CSV-Datei einlesen (mit pandas)
def read_csv_pandas(filename):
    df = pd.read_csv(filename, header=None)  # Liest CSV ohne Kopfzeile
    df = df.dropna(axis=1)  # Entfernt leere Spalten (falls welche existieren)
    return df.values.flatten().astype(float)  # Konvertiert in flaches NumPy-Array (float, um Fehler zu vermeiden)

# Daten einlesen
data = read_csv_pandas("outcomes5.csv")
print(str(len(data)))

# Plot erstellen
plt.figure(figsize=(8, 5))
plt.plot(data, marker="o", linestyle="-", color="b", linewidth=1, markersize=1)

# Styling für einen coolen Look
plt.xlabel("n/10", fontsize=12)
plt.ylabel("p(FLOPS)", fontsize=12)
plt.xscale('log')
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.show()

