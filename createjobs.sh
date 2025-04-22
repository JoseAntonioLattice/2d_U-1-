#!/bin/bash


filename=2dU1_L100_
jobname=L100_
#touch $filename$1.slurm

cat <<EOF >> jobs/$filename$1.slurm
#!/bin/bash
#SBATCH --job-name=logs/$jobname$1 			# Nombre del trabajo
#SBATCH --output=logs/$jobname$1.log   		# Archivo de registro de salida
#SBATCH --error=logs/$jobname$1.err    		# Archivo de registro de errores
#SBATCH --partition=QuantPhysMC   	# Nombre de la partición o cola de trabajos
#SBATCH --nodes=1     				# Número de nodos a utilizar (puedes cambiarlo)
#SBATCH --ntasks-per-node=1   		# Número de tareas por nodo (1 para ejecución serial)
#SBATCH --cpus-per-task=1 		# Número de CPUs por tarea (puedes cambiarlo)
#SBATCH --mem=4G      			# Memoria RAM necesaria (puedes cambiarlo)

module load lamod/gcc/12.2 
cd /home/icn/joseantog/2d_U-1-
# Comando para ejecutar tu programa
make run

EOF
