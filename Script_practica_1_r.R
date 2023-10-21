#############################################################################
#
# PRACTICA 3
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
## ENTREGA EL 22 OCTUBRE 23:59
## Entrega de este script con los códigos más las imágenes y las respuestas a las preguntas
## Adjuntar en la entrega el PDF final y el archivo con los genes
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")
## nota: no ha habido problema

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))
he
# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data)
##nota: sirve para medir dimensiones
head (data)
tail (data)

# Hacemos un primer histograma para explorar los datos
hist(data)
hist (data, col="green")

# Transformamos los datos con un logaritmo 
histogram= log(data)
#CAMBIO DE NOMBRE A LA VARIABLE
hist (histogram)
hist (histogram, col= "purple")
##nota: no dejamos espacio entre log y (data)
## guardamos la transformacion logaritmica en forma de variable

# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
##respuesta: La transformación logaritmica se emplea para que los datos sean más fáciles 
de interpretar (campana de Gauss) y para hacer imágenes que puedan ser publicadas en artículos científicos

# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
boxplot(histogram)
boxplot(histogram, col=c("blue","blue","blue","orange","orange","orange"),main="boxplotdata",las=2)
##nota: le ponemos colores diferentes 
##nombrado como "boxplotdata"
##las=2 --> permite poner los nombres identificadores en vertical 

# ¿Qué es un boxplot?
##respuesta: un boxplot es una representación gráfica de una serie de datos que se disponen en barras
alineadas según la media, la mediana y los cuartiles. Este tipo de gráfico se caracteriza por presentar bigotes que señalan 
los datos máximos y mínimos

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación
# de los valores de expresión. ¿Es correcta la separación?
hc=hclust(as.dist(1-cor(histogram)))
plot(hc, main="clustering")
##nota: una coreccta clasificación es si tenemos 3 grupos de WT y otros 3 grupos de mutados
##respuesta: sí

#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
wt <- data [,1:3]
ko <- data [,4:6]
class(wt)
##nota: indica el tipo de datos que tenemos(una matriz o tabla)
head(wt)
##respuesta: hemos generado datos en forma de tabla o matriz

# Calcula las medias de las muestras para cada condición. Usa apply (podemos calcular donde queramos)
## creamos una variable 
wt.mean = apply(wt, 1, mean)
ko.mean = apply(ko, 1, mean)
##nota: queremos ver si la expresión de los ko son difentes a los wt

# ¿Cuál es la media más alta?
max(wt.mean)
##media más alta de wt 
max(ko.mean)
## respuesta: media más alta es la del grupo de ko

# Ahora hacemos un scatter plot (gráfico de dispersión)
plot(ko.mean ~ wt.mean)
## nota: compara medias en el eje x y en el eje y
plot(ko.mean ~ wt.mean, xlab= "WT", ylab="KO", main="scatterplot")

# Añadir una línea diagonal con abline
## poner una linea de regresión. lo hacemos con abline, se hace cuando el plot está abierto
abline (0, 1, col= "red")
## podemos poner una linea horizontal, ponemos una h (horizontal)para una linea horizontal, o vertical con v
abline (h=2, col="blue")
abline (v=5, col="green")

# ¿Eres capaz de añadirle un grid?
grid ()
## añadimos cuadricula

# Calculamos la diferencia entre las medias de las condiciones
diff.mean = wt.mean - ko.mean
## hemos obtenido una tabla (porque guardanos todos los datos)

# Hacemos un histograma de las diferencias de medias
hist( diff.mean)
hist (diff.mean, col="pink", xlab="diferencia de medias")

# Calculamos la significancia estadística con un t-test.
##nota: hemos empleado un bucle tipo for
pvalue=NULL
tstat=NULL
for (i in 1 : nrow(data)) { #Para cada gen 
	x = wt[i,] #gene wt nº i
	y = ko[i,] #gene ko nº i
##hacemos el test 
	t = t.test(x,y)
##añadimos el p-value a la lita 
	pvalue[i] = t$p.value
##añadimos las estadísticas a la lista
	tstat[i]= t$statistic
}
head(pvalue)

length (pvalue)
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
## respuesta:para hacer el análisis utilizamos los datos obtenidos sin transformar en
logaritmo para evitar modificar los datos y alterar los resultados. 

# ¿Cuántas valores tiene cada muestra?
##tenemos 3 muestras para cada dos condiciones

# Ahora comprobamos que hemos hecho TODOS los cálculos
##comprobado: sí

# Hacemos un histograma de los p-values.
hist (pvalue)

# ¿Qué pasa si le ponemos con una transformación de -log10?
hist (-log10(pvalue), col="blue")
##respuesta: al poner los datos con una transformación logaritmica de -log10, los datos se muestran más ordenados y homogéneos en el histograma, lo que 
facilita su interpretación

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(diff.mean,-log10(pvalue), main=" Volcano - pvalue", col="orange", xlab="diferencia de medias en contra del log")
## nota: pongo el pvalue que es lo que quiero representar. los pvalue significativos, 
se encuentran en la cola del volcano

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
diff.mean_cutoff = 2
pvalue_cutoff = 0.01
abline (v = diff.mean_cutoff, col="blue", lwd = 3)
## ponemos los pvalue en la linea horizontal 
abline(h= -log10(pvalue_cutoff), col="green", lwd= 3)
##marcamos a partir de donde hay significacion

# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean = abs (diff.mean)>= diff.mean_cutoff
dim(data[filter_by_diff.mean])
## nota: convierto en valores absolutos las diferencias de medias por encima de 2

# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data[filter_by_pvalue,])
##extraemos todos los pvalues menores o iguales. pongo el dim para buscar cuantos genes 
sobrepasan el filtro

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?
filter_combined = filter_by_diff.mean & filter_by_pvalue
filtered = data [filter_combined,]
dim (filtered)
head (filtered)
## combino ambos identificadores y me quedo con los valores que me den los filtros
extraemos todo

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean,-log10(pvalue), main ="Volcano con valores significativos")
points(diff.mean[filter_combined], -log10(pvalue[filter_combined]),col= "red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés?
##respuesta: al hacer la diferencia de medias calculando wt-ko , los valores del ko son mayores, por lo que 
la diferencia de medias es negativa y los sobrexpresados aparecen en el cuadrante negativo
plot(diff.mean, -log10(pvalue), main="volcano 3")
plot(diff.mean[filter_combined & diff.mean < 0], 
- log10(pvalue[filter_combined & diff.mean <0], col="red")
points (diff.mean[filter_combined & diff.mean >0], 
-log10(pvalue[filter_combined & diff.mean > 0]) col="blue")

# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7, labRow=FALSE)
##agrupa los wt y ko, y los coloca en la mitad en sobreexpresados y reprimidos
heatmap(filtered)


# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)

# Hacemos nuestro heatmap
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7, 
col=rev(redblue(256)), scale="row")

# Lo guardamos en un archivo PDF


# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table(filtered, "genes_practica_1R.txt",sep="\t", quote=FALSE)
