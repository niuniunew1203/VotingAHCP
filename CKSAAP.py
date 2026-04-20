import csv


def CKSAAP(seq):
    from itertools import product

    allowed_chars = 'ACDEFGHIKLMNPQRSTVWY'

    total_data = []
    for i in seq:
        if i not in allowed_chars:
            print("Invalid string")
            exit()

    object_data = []

    for k in range(4):
        pairs = [''.join(pair) for pair in product(allowed_chars, repeat=2)]
        pairs = {pair: 0 for pair in pairs}
        for i in range(len(seq)):
            for j in range(i + k + 1, len(seq)):
                if j - i == k + 1:
                    pair = seq[i] + seq[j]
                    if set(pair).issubset(set(allowed_chars)):
                        if pair in pairs:
                            pairs[pair] += 1
                        else:
                            pairs[pair] = 1
        if k == 0:
            for l in pairs:
                pairs[l] = pairs[l] / (len(seq) - 1)
        elif k == 1:
            for l in pairs:
                pairs[l] = pairs[l] / (len(seq) - 2)
        elif k == 2:
            for l in pairs:
                pairs[l] = pairs[l] / (len(seq) - 3)
        elif k == 3:
            for l in pairs:
                pairs[l] = pairs[l] / (len(seq) - 4)

        object_data.extend([pairs[pair] for _ in range(1) for pair in pairs])
        # print(len(object_data))

    # print(object_data)

    total_data.append([seq] + object_data)

    # print(total_data)
    # print(len(total_data[0]))

    headers = ['Seq'] + [f'k{k}_{pair}' for k in range(4) for pair in pairs]

    filename = "CKSAAP.csv"
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerow(total_data[0])  # Write the row containing all the data

    print(f"CSV file '{filename}' created successfully!")


seq = 'VSFAIKWEYVLLL'
CKSAAP(seq)
