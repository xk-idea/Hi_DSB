import os

def count_labels_in_file(file_path):

    count_0 = 0
    count_1 = 0
    try:
        with open(file_path, 'r') as f:
            for line in f:

                labels = line.strip().split()
                for label in labels:
                    if label == '0':
                        count_0 += 1
                    elif label == '1':
                        count_1 += 1
        return count_0, count_1
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return 0, 0
    except Exception as e:
        print(f"Error occurred while reading {file_path} : {e}")
        return 0, 0

def main():
    directory = "../preprocess/Label"
    
    total_0 = 0
    total_1 = 0

    for chrom in range(1, 23):
        filename = f"chr{chrom}_label.txt"
        file_path = os.path.join(directory, filename)
        count_0, count_1 = count_labels_in_file(file_path)
        print(f"{filename} - 0: {count_0}, 1: {count_1}")
        total_0 += count_0
        total_1 += count_1

    print(f"Count0: {total_0}")
    print(f"Count1: {total_1}")

if __name__ == "__main__":
    main()