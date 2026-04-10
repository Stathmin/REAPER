#!/bin/bash
# Initialize Probename project with genome sizes and comparative analyses
# Follows holy principles: modular architecture, configuration compliance, automated project management

echo "🔬 Initializing Probename project with genome size-based comparative analysis..."

# Check if project manager exists
if [ ! -f "project_manager.py" ]; then
    echo "❌ project_manager.py not found"
    exit 1
fi

# Check if data files exist
if [ ! -f "data/probename_reads/testicum_R1.fq" ] || [ ! -f "data/probename_reads/testicum_R2.fq" ]; then
    echo "❌ testicum data not found: data/probename_reads/testicum_R1.fq"
    exit 1
fi

if [ ! -f "data/probename_reads/examicum_R1.fq" ] || [ ! -f "data/probename_reads/examicum_R2.fq" ]; then
    echo "❌ examicum data not found: data/probename_reads/examicum_R1.fq"
    exit 1
fi

if [ ! -f "data/probename_reads/studicum_R1.fq" ] || [ ! -f "data/probename_reads/studicum_R2.fq" ]; then
    echo "❌ studicum data not found: data/probename_reads/studicum_R1.fq"
    exit 1
fi

# 🏛️ HOLY PRINCIPLE: Configuration Compliance - Clean slate approach
echo "🧹 Ensuring clean slate for project initialization..."
if [ -d "projects/probename_project" ]; then
    echo "📁 Removing existing probename_project to prevent interference..."
    rm -rf projects/probename_project
    echo "✅ Clean slate prepared"
fi

# Create project
echo "📁 Creating probename_project..."
python3 project_manager.py create-project \
    --project-id probename_project \
    --taxonomy "probename" \
    --description "Probename family repeat analysis with genome size-based comparative analysis" \
    --ncbi-repeats "data/ncbi_repeats.fa" \
    --comparative-species "testicum" "examicum" "studicum"

if [ $? -ne 0 ]; then
    echo "❌ Failed to create project"
    exit 1
fi

# Add samples with genome sizes
echo "🧬 Adding testicum sample (genome_size: 1.0)..."
python3 project_manager.py add-sample \
    --project-id probename_project \
    --sample-id testicum \
    --taxonomy "probename" \
    --r1-path "data/probename_reads/testicum_R1.fq" \
    --r2-path "data/probename_reads/testicum_R2.fq" \
    --genome-size 1.0

if [ $? -ne 0 ]; then
    echo "❌ Failed to add testicum sample"
    exit 1
fi

echo "🧬 Adding examicum sample (genome_size: 1.2)..."
python3 project_manager.py add-sample \
    --project-id probename_project \
    --sample-id examicum \
    --taxonomy "probename" \
    --r1-path "data/probename_reads/examicum_R1.fq" \
    --r2-path "data/probename_reads/examicum_R2.fq" \
    --genome-size 1.2

if [ $? -ne 0 ]; then
    echo "❌ Failed to add examicum sample"
    exit 1
fi

echo "🧬 Adding studicum sample (genome_size: 1.5)..."
python3 project_manager.py add-sample \
    --project-id probename_project \
    --sample-id studicum \
    --taxonomy "probename" \
    --r1-path "data/probename_reads/studicum_R1.fq" \
    --r2-path "data/probename_reads/studicum_R2.fq" \
    --genome-size 1.5

if [ $? -ne 0 ]; then
    echo "❌ Failed to add studicum sample"
    exit 1
fi

# Verify project structure
echo "🔍 Verifying project structure..."
python3 project_manager.py validate --project-id probename_project

if [ $? -ne 0 ]; then
    echo "❌ Project validation failed"
    exit 1
fi

# 🏛️ HOLY PRINCIPLE: Update existing samples with unique prefixes
echo "🔤 Updating samples with unique 4-letter prefixes..."
python3 project_manager.py update-prefixes --project-id probename_project

if [ $? -ne 0 ]; then
    echo "❌ Failed to update prefixes"
    exit 1
fi

echo "✅ Unique prefixes assigned to all samples"

# Add comparative analysis automatically
echo "🔬 Adding comparative analysis: testicum vs examicum..."
python3 project_manager.py add-comparative \
    --project-id probename_project \
    --analysis-id comp_testicum_examicum \
    --samples testicum examicum \
    --analysis-description "Testicum vs Examicum comparison"

if [ $? -ne 0 ]; then
    echo "❌ Failed to add comparative analysis"
    exit 1
fi

echo "✅ Probename project initialized successfully!"
echo ""
echo "📊 Project Summary:"
echo "  - Project: probename_project"
echo "  - Samples: testicum (genome_size: 1.0), examicum (genome_size: 1.2), studicum (genome_size: 1.5)"
echo "  - Comparative Analysis: comp_testicum_examicum (testicum vs examicum)"
echo ""
echo "📋 Available comparative analyses:"
python3 project_manager.py list-comparatives --project-id probename_project
echo ""
echo "🚀 Ready to run pipeline with:"
echo "  snakemake -s Snakefile_modular --configfile projects/global_config.yaml --cores 4"
echo ""
echo "🎯 Pipeline targets:"
echo "  - Individual assemblies: testicum, examicum, studicum"
echo "  - Comparative analysis: comp_testicum_examicum"
echo ""
echo "📊 To run all analyses:"
echo "  snakemake -s Snakefile_modular --configfile projects/global_config.yaml --cores 4 --target all"
