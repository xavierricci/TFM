# TFM
Creación de un modelo de predicción de riesgos de incendios forestales usando una red neuronal convolucional sobre datos históricos de meteorología de California.


Codigos creados y resultados

Este conjunto de códigos es parte del TFM de Xavier Ricci y repasa los procesos seguidos para la generación de una red neuronal convolucional que, a partir del histórico de los incendios sucedidos en california entre 2000 y 2019; y usando datos meteorologicos y topograficos de la fecha y lugar de cada uno de estos incendios, es capaz de predecir las zonas con mayor riesgo de incendios teniendo en cuenta la situacion climatologica.  

Este sistema de predicción se complementa con un generador de heatmans que a partir de unas coordenadas y datos metereologicos es capaz de generar un heatmap a 500m de resolucion sobre las probabilidades de incendio de una zona especificada.  

El sistema ha sido testado con datos de fuegos de 2019 comparados con datos de puntos random con vegetacion de California y ha acertado en un 95% que punto fue un incendio y que punto no lo fué.
