import random

def generate_structure(structure_length=400, open_region_length=50):
    region = []
    while len(region) < open_region_length:
        region_length = int(random.random()*10) + 5
        region_start = int(random.random()*structure_length)
        region_end = region_start + region_length
        if region_start - 5 in region or region_end + 5 in region or region_start < 0 or region_end > structure_length:
            continue
        for i in range(region_start,region_end+1):
            region.append(i)
    return sorted(region)


def generate_read(structure, open_regions, read_length=400, mod_ratio=0.10):
    """
    生成单个read的序列信息
    """
    read = [0] * read_length
    for i in open_regions:
        if random.random() < mod_ratio:
            read[i] = 1
    return read



def generate_reads(num_reads, structure, open_regions, read_length=400, mod_ratio=0.10):
    """
    生成多条reads的序列信息
    """
    return [generate_read(structure, open_regions, read_length, mod_ratio) for _ in range(num_reads)]





structure1_open_regions = generate_structure()
structure2_open_regions = generate_structure()

# 生成结构1和结构2的reads
structure1_reads = generate_reads(500, "Structure1", structure1_open_regions)
structure2_reads = generate_reads(500, "Structure2", structure2_open_regions)

print(f"#{structure1_open_regions}")
print(f"#{structure2_open_regions}")
positive_ratio= 0.02
for n in range(500):
    index = f'structure1_read_{n}'
    for i in range(400):
        if random.random() < positive_ratio:
            structure1_reads[n][i] = 1
        print(f'{index}\t{i}\t{structure1_reads[n][i]}')
        
for n in range(500):
    index = f'structure2_read_{n}'
    for i in range(400):
        if random.random() < positive_ratio:
            structure1_reads[n][i] = 1
        print(f'{index}\t{i}\t{structure2_reads[n][i]}')
