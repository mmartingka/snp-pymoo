import os
import re
import pandas as pd
import argparse

def procesar_archivos_snps(directorio_base, salida="resultado.csv"):
    # Crear una lista para almacenar los datos extraídos
    datos = []

    # Expresión regular para capturar los tres últimos números en el nombre del archivo (pueden ser negativos)
    patron = re.compile(r".*-(\d+)-(\d+)-(\d+)-PF\.csv$")

    # Recorrer todas las carpetas y subcarpetas
    for carpeta_inicial in os.listdir(directorio_base):
        ruta_carpeta_inicial = os.path.join(directorio_base, carpeta_inicial)
            
        # Verifica que sea un directorio
        if not os.path.isdir(ruta_carpeta_inicial):
            continue

        for subcarpeta_hip in os.listdir(ruta_carpeta_inicial):
            # Solo procesar las carpetas que comienzan con "hip"
            if not subcarpeta_hip.startswith("hip"):
                continue
            
            ruta_subcarpeta = os.path.join(ruta_carpeta_inicial, subcarpeta_hip)
            
            # Verifica que sea un directorio
            if not os.path.isdir(ruta_subcarpeta):
                continue
            
            for archivo in os.listdir(ruta_subcarpeta):
                # Verificar que el archivo tenga extensión .csv
                if archivo.endswith(".csv"):
                    ruta_archivo = os.path.join(ruta_subcarpeta, archivo)
                    print

                    print(ruta_archivo)
                    # Buscar los tres últimos hiperparámetros en el nombre del archivo
                    coincidencia = patron.search(archivo)
                    print(coincidencia)
                    if coincidencia:
                        # Convertir los valores a enteros
                        l_mut, f_mut, prob_cross = coincidencia.groups()
                    else:
                        print(f"Nombre de archivo no válido o sin hiperparámetros reconocibles: {archivo}")
                        continue

                    if '--' in archivo:
                        l_mut = '-'+l_mut

                    # Leer el archivo y calcular la media
                    df = pd.read_csv(ruta_archivo, header=None)
                    valor = df[0][0]
                    
                    # Guardar la información en la lista de datos
                    datos.append({
                        "Directorio": ruta_carpeta_inicial,
                        "Subcarpeta": subcarpeta_hip,
                        "L_mut": l_mut,
                        "F_mut": f_mut,
                        "Prob_cross": prob_cross,
                        "Media": float(valor)
                    })
    
    # Convertir los datos en un DataFrame
    df_resultado = pd.DataFrame(datos)

    # Agrupar por hiperparámetros y calcular la media general de los valores para los mismos
    df_agrupado = df_resultado.groupby(["Directorio", "L_mut", "F_mut", "Prob_cross"]).agg({
        "Media": "mean"
    }).reset_index()

    # Ordenar los valores por proximidad a 1.000 en la columna "Media"
    df_ordenado = df_agrupado.sort_values(by="Media", key=lambda x: abs(x - 1.0)).head(10)

    # Guardar el resultado en un archivo CSV
    df_ordenado.to_csv(salida, index=False)
    print(f"Archivo de resultados guardado como '{salida}'")

# Ejecutar el script si se llama desde la línea de comandos
if __name__ == "__main__":
    # Parsear argumentos de línea de comandos
    parser = argparse.ArgumentParser(description="Procesar archivos de SNPs en directorios y subdirectorios.")
    parser.add_argument("directorio", type=str, help="Directorio base que contiene las carpetas de SNPs")
    parser.add_argument("--salida", type=str, default="resultado.csv", help="Nombre del archivo de salida (por defecto: resultado.csv)")

    args = parser.parse_args()

    # Llamar a la función con el directorio proporcionado
    procesar_archivos_snps(args.directorio, args.salida)