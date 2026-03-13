# Leemos los datos de consumo de glucosa, produccion de acetato, lactato, formiato y crecimiento en OD
# Creamos una funcion de lectura de datos de consumo/produccion de metabolitos, el formato de la tabla es:

# #: es el indice
# Condicion: nombre de la condicion, puede ser control, WT, 
# Muestra: nombre de la muestra, puede ser t:0 control, ΔlppB, ΔlppA
# time[h]: tiempo de la medicion
# OD 600: valor de la OD en la medicion
# Glucosa: valor de la concentracion de la glucosa en la medicion
# Formic acid: valor de la concentracion de formiato en la medicion
# Lactic acid: valor de la concentracion de lactato en la medicion
# Acetic acid: valor de la concentracion de acetato en la medicion

def read_metabolite_consumption_production(file_path, metabolite_columns=["Glucosa", "Formic acid", "Lactic acid", "Acetic acid"], od_column="OD600", time_column="time[h]", condition_column="Condition"):
    # Leemos los datos de abundancia de metabolitos, devolviendo un dataframa. Los valores vacios se reemplazan por 0.
    data_df=pd.read_excel(file_path,sheet_name="datos")
    data_df=data_df.fillna(0)

    return(data_df)

# Dado que se leen abundancias de metabolitos, vamos a crear una funcion para transformar estas abundancias en tasas de consumo/produccion, para esto necesitamos la OD y el tiempo de cada medicion, ademas de un factor de conversion entre OD y biomasa (gDW/L). El volumen del cultivo (L) será 1.
# Como entrada tendra un dataframe generado por read_metabolite_consumption_production ya añadira columnas para el calculo de las tasas, devolviendo un dataframe con las mismas columnas pero con las tasas de consumo/produccion en lugar de las abundancias. La formula para el calculo de las tasas sera: (abundancia final - abundancia inicial) / (biomasa * tiempo), asumiendo que la abundancia inicial es 0.02.
# Referencia: Raghunathan, A., Reed, J., Shin, S., Palsson, B., & Daefler, S. (2009). Constraint-based analysis of metabolic capacity of Salmonella typhimurium during host-pathogen interaction. BMC systems biology, 3(1), 38.
def abundance_to_rate(OD_ini, OD_fin, time_ini,time_fin, abundance_ini, abundance_fin, conversion_factor=0.7396):
    # Calculamos la tasa de consumo/produccion de un metabolito a partir de su abundancia inicial y final, el tiempo entre las mediciones y la OD inicial y final. La formula para el calculo de las tasas sera: (abundancia final - abundancia inicial) / (biomasa * tiempo), asumiendo que la abundancia inicial es 0.02.
    biomass_ini=OD_ini*conversion_factor
    biomass_fin=OD_fin*conversion_factor
    biomass_avg=(biomass_ini+biomass_fin)/2
    delta_time=time_fin-time_ini
    if abundance_ini>abundance_fin:
        delta_abundance=abundance_ini-abundance_fin
        delta_abundance=-delta_abundance # Si la abundancia final es menor que la inicial, es un consumo, por lo que la tasa debe ser negativa
    else:
        delta_abundance=abundance_fin-abundance_ini
    if delta_time==0 or biomass_avg==0:
        rate=0
    else:
        rate=(delta_abundance)/(biomass_avg*delta_time)
    return(rate)

def transform_abundance_to_rates(abundance_df, condition_column="Condicion", od_column="OD600",metabolite_columns=["Glucosa", "Formic acid", "Lactic acid", "Acetic acid"],time_column="time[h]", conversion_factor=0.7396):
    # Transformamos las abundancias de metabolitos en tasas de consumo/produccion, devolviendo un dataframe con las mismas columnas pero con las tasas de consumo/produccion en lugar de las abundancias. La formula para el calculo de las tasas sera: (abundancia final - abundancia inicial) / ((biomasa media)* (tiempo final tiempo inicial)), asumiendo que la abundancia inicial es el valor del tiempo anterior siendo el primer valor 0.02. El volumen del cultivo (L) será 1.
    # Agrupaos los datos por condicion para aplicar la transformacion de abundancia a tasas por separado para cada condicion
    abundance_df=abundance_df.copy()
    grouped = abundance_df.groupby(condition_column)
    # Cogemos las filas del grupo control para obtener los valores iniciales de tiempo, OD y abundancia de metabolitos, asumiendo que el control es la condicion que tiene los valores iniciales de tiempo, OD y abundancia de metabolitos. Si no hay un grupo llamado control, se tomara el primer grupo como control.
    control=grouped.get_group("control")
    time_0=control[time_column].min()
    OD_0=control[od_column].min()
    metabolites_t0={}
    for metabolite in metabolite_columns:
        metabolites_t0[metabolite]=control[metabolite].min()
    for name, group in grouped:
        group = group.sort_values(by=time_column)
        for metabolite in metabolite_columns:
            rates = []
            for i in range(len(group)):
                if i == 0:
                    # Si es el primer tiempo, poner los valores de time_0 que tiene los valores iniciales de abundancia y OD
                    rate = abundance_to_rate(OD_ini=OD_0, OD_fin=group.iloc[i][od_column], time_ini=time_0, time_fin=group.iloc[i][time_column], abundance_ini=metabolites_t0[metabolite], abundance_fin=group.iloc[i][metabolite], conversion_factor=conversion_factor)
                else:
                    rate = abundance_to_rate(OD_ini=group.iloc[i-1][od_column], OD_fin=group.iloc[i][od_column], time_ini=group.iloc[i-1][time_column], time_fin=group.iloc[i][time_column], abundance_ini=group.iloc[i-1][metabolite], abundance_fin=group.iloc[i][metabolite], conversion_factor=conversion_factor)
                rates.append(rate)
            # Creamos una nueva columna en el dataframe original para almacenar las tasas de consumo/produccion de cada metabolito, con el nombre del metabolito seguido de "_rate"
            abundance_df.loc[group.index, metabolite + "_rate"] = rates
    
    return abundance_df.reset_index(drop=True)

def calculate_growth_rate(abundance_df, condition_column="Condicion", od_column="OD600", time_column="time[h]", conversion_factor=0.7396):
    # Calculamos la tasa de crecimiento a partir de la OD y el tiempo, devolviendo un dataframe con las mismas columnas pero con una nueva columna "growth_rate" que contiene la tasa de crecimiento calculada. La formula para el calculo de la tasa de crecimiento sera: (OD_final - OD_inicial) / (tiempo_final - tiempo_inicial) / conversion_factor, asumiendo que la OD inicial es el valor del tiempo anterior siendo el primer valor 0.02.
    abundance_df=abundance_df.copy()
    grouped = abundance_df.groupby(condition_column)
    control=grouped.get_group("control")
    time_0=control[time_column].min()
    OD_0=control[od_column].min()
    for name, group in grouped:
        group = group.sort_values(by=time_column)
        growth_rates = []
        for i in range(len(group)):
            if i == 0:
                growth_rate = (group.iloc[i][od_column] - OD_0) / (group.iloc[i][time_column] - time_0) / conversion_factor
            else:
                growth_rate = (group.iloc[i][od_column] - group.iloc[i-1][od_column]) / (group.iloc[i][time_column] - group.iloc[i-1][time_column]) / conversion_factor
            growth_rates.append(growth_rate)
        abundance_df.loc[group.index, "growth_rate"] = growth_rates
    
    return abundance_df.reset_index(drop=True)

data_dict=read_metabolite_consumption_production("proteomics/HPLC_results-processed_Juanjo_2.xlsx")
data_rates_df=transform_abundance_to_rates(data_dict)
data_rates_df=calculate_growth_rate(data_rates_df)
