home_dir = config['home_dir']

tools = {
    'qctool':f"{home_dir}software/QCTool/qctool/build/release/apps/qctool_v2.2.2",
    'snptest': f"{home_dir}software/snptest/snptest_v2.5.6",
    'inthinnerator': f"{home_dir}software/QCTool/qctool/build/release/apps/inthinnerator_v2.2.2",

    'picard': f"java -Xmx20G -Xms10G -jar {home_dir}conda/skylake/envs/lcwgs/share/picard-slim-2.27.4-0/picard.jar",
    'picard_plus': f"java -Xmx40G -Xms20G -jar {home_dir}conda/skylake/envs/lcwgs/share/picard-slim-2.27.4-0/picard.jar",
    'picard_pplus': f"java -Xmx60G -Xms30G -jar {home_dir}conda/skylake/envs/lcwgs/share/picard-slim-2.27.4-0/picard.jar",
    'picard_ppplus': f"java -Xmx100G -Xms40G -jar {home_dir}conda/skylake/envs/lcwgs/share/picard-slim-2.27.4-0/picard.jar",
    'gatk': "gatk --java-options -Xmx8G",

    'quilt_hla':f"{home_dir}software/QUILT/QUILT_HLA.R",
    'quilt_hla_prep': f"{home_dir}software/QUILT/QUILT_HLA_prepare_reference.R",
    
    'hla-la': "/well/ansari/users/gjx698/software/HLA-LA/src/HLA-LA.pl",
    'impute2': f"{home_dir}software/impute2/impute2",
    'beagle': f"java -Xmx60G -Xms30G -jar {home_dir}software/BEAGLE/beagle.29Oct24.c8e.jar",
    'coverotron': f'{home_dir}software/dir/build/apps/coverotron',
    'classify-kmers': f'{home_dir}software/dir/build/apps/classify-kmers'
}
