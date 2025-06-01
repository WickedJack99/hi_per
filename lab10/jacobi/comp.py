import csv
import math

def compare_csv(file1, file2):
    with open(file1, newline='') as f1, open(file2, newline='') as f2:
        reader1 = list(csv.reader(f1))
        reader2 = list(csv.reader(f2))

        if len(reader1) != len(reader2):
            print(f"Unterschiedliche Zeilenzahl: {len(reader1)} vs {len(reader2)}")
            return

        for row_idx, (row1, row2) in enumerate(zip(reader1, reader2), start=1):
            if len(row1) != len(row2):
                print(f"Zeile {row_idx}: unterschiedliche Spaltenanzahl")
                continue
            for col_idx, (cell1, cell2) in enumerate(zip(row1, row2), start=1):
                if not math.isclose(float(cell1), float(cell2), abs_tol=1e-7):
                    print(f"Unterschied in Zeile {row_idx}, Spalte {col_idx}: {cell1} ≠ {cell2} diff: {abs(float(cell1) - float(cell2))}")

if __name__ == "__main__":
    compare_csv("result_parallel.asc", "result_serial.asc")

