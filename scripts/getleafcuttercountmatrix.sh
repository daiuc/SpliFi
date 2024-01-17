#!/bin/bash

echo "## $(date): Getting Leafcutter count matrices for all tissues"

tissues='Adipose-Subcutaneous Adipose-Visceral_Omentum_ AdrenalGland Artery-Aorta Artery-Coronary Artery-Tibial Bladder Brain-Amygdala Brain-Anteriorcingulatecortex_BA24_ Brain-Caudate_basalganglia_ Brain-CerebellarHemisphere Brain-Cerebellum Brain-Cortex Brain-FrontalCortex_BA9_ Brain-Hippocampus Brain-Hypothalamus Brain-Nucleusaccumbens_basalganglia_ Brain-Putamen_basalganglia_ Brain-Spinalcord_cervicalc-1_ Brain-Substantianigra Breast-MammaryTissue Cells-Culturedfibroblasts Cells-EBV-transformedlymphocytes Colon-Sigmoid Colon-Transverse Esophagus-GastroesophagealJunction Esophagus-Mucosa Esophagus-Muscularis Heart-AtrialAppendage Heart-LeftVentricle Kidney-Medulla Liver Lung MinorSalivaryGland Muscle-Skeletal Nerve-Tibial Ovary Pancreas Pituitary Prostate Skin-NotSunExposed_Suprapubic_ Skin-SunExposed_Lowerleg_ SmallIntestine-TerminalIleum Spleen Stomach Testis Thyroid Uterus Vagina WholeBlood'

for t in ${tissues[@]}; do 

    echo Rscript getleafcuttercountmatrix.R "$t"
    Rscript getleafcuttercountmatrix.R "$t" &

    if (( $(jobs | wc -l) > 5 )); then
        echo "wait for previous jobs to finish"
        wait    # Wait for any of the background processes to finish
        echo $t done at $(date)
    fi
    
done

echo "## $(date): all done."

