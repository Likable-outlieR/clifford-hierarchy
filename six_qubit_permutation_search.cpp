#include <iostream>
using namespace std;

int dig_sum(int k) {
    // computes the sum of the binary digits of a number k
    return (k ? (dig_sum(k>>1) + (k&1)) : 0);
}

bool is_clifford_perm(int n, int* sigma) {
    // given an array sigma (assumed to be of length 2 to the n), check whether its entries are a clifford permutation
    int N = 1<<n;

    int* shifted_sigma = new int[N];

    for(int j=0; j<N; j++) {
        shifted_sigma[j] = sigma[j]^sigma[0];
    }

    bool ans = true;

    for(int j=0; j<N; j++) { //what follows uses a trick from folklore: for any positive integer j, j&(j-1) is j with its last 1 removed
        if(shifted_sigma[j] != (shifted_sigma[j&(j-1)]^shifted_sigma[j-(j&(j-1))])) {
            ans = false;
        }
    }

    return ans;
}

bool yields_x(int n, int x, int* sigma, int* sigma_inv) {
    // given an array sigma, assumed to be a permutation whose inverse is sigma_inv,
    // checks whether the result of conjugating [the pauli-x operator indexed by the binary digits of x] by sigma
    // is a pauli-x operator
    int N = 1<<n;

    int c0 = sigma[sigma_inv[0]^x];

    bool ans = true;

    for(int j=1; j<N; j++) {
        if(sigma[sigma_inv[j]^x] != (c0^j)) { // sigma[sigma_inv[j]^x] gives the new conjugated permutation
            ans = false;
        }
    }

    return ans;
}

bool yields_z(int n, int z, int* sigma, int* sigma_inv) {
    // given an array sigma, assumed to be a permutation whose inverse is sigma_inv,
    // checks whether the result of conjugating [the pauli-z operator indexed by the binary digits of z] by sigma
    // is a pauli-z operator
    int N = 1<<n;

    int* conj_z = new int[N];

    for(int j=0; j<N; j++) {
        conj_z[j] = (dig_sum(sigma_inv[j]&z)&1);
    }
    // the result of conjugating [the pauli-z operator indexed by the binary digits of z] by sigma
    // is diagonal with 1's and -1's. the 1's and -1's are given by the 0's and 1's of conj_z.
    // so now we check whether conj_z defines a pauli-z gate

    for(int j=1; j<N; j++) {
        conj_z[j] ^= conj_z[0];
    }

    conj_z[0] = 0;

    bool ans = true;

    for(int j=0; j<N; j++) { //this is the same trick as above
        if(conj_z[j] != (conj_z[j&(j-1)]^conj_z[j-(j&(j-1))])) {
            ans = false;
        }
    }

    return ans;
}

void print_bin(int n, int x) {
// prints a n-digit number x in binary. leading zeroes are allowed
    for(int i=1; i<=n; i++) {
        if(x&(1<<(n-i))) {
            cout << "1";
        } else {
            cout << "0";
        }
    }

    return;
}

int main() {
    // the controls and targets of all possible toffoli gates in a six-qubit staircase form permutation
    int control_one[20] = {1,1,1,2,1,1,2,1,2,3,1,1,2,1,2,3,1,2,3,4};
    int control_two[20] = {2,2,3,3,2,3,3,4,4,4,2,3,3,4,4,4,5,5,5,5};
    int target[20]      = {3,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6};

    cout << "starting search" << '\n' << '\n';
    //int count = 0;

    for(int i=0; i<(1<<20); i++) {  // we search over all 2 to the 20 staircase form permutations pi on 6 qubits.
                                    // the binary digits of i tell which toffoli gates to include
        int pi[64];

        for(int j=0; j<64; j++) {
            pi[j] = j;
        }

        for(int k=0; k<20; k++) {
            if(i&(1<<k)) {
                int c1 = 6-control_one[k]; //the indexing is done like this because we number the qubits left to right
                int c2 = 6-control_two[k];
                int t  = 6-target[k];

                for(int j=0; j<64; j++) {
                    if((pi[j]&(1<<c1))&&(pi[j]&(1<<c2))) {
                        pi[j] = pi[j]^(1<<t);
                    }
                }
            }
        }

        int pi_inv[64];

        for(int j=0; j<64; j++) {
            pi_inv[pi[j]] = j;
        }

        bool is_in_c3 = true;

        for(int m=0; m<6; m++) {
            int tau[64];

            for(int j=0; j<64; j++) {
                tau[j] = pi[pi_inv[j]^(1<<m)];
            }

            if(!(is_clifford_perm(6, tau))) {
                is_in_c3 = false;
            }
        }

        bool control_used[7];
        bool target_used[7];

        for(int p=0; p<=6; p++) {
            control_used[p] = false;
            target_used[p] = false;
        }

        for(int k=0; k<20; k++) {
            if(i&(1<<k)) {
                control_used[control_one[k]] = true;
                control_used[control_two[k]] = true;
                target_used[target[k]] = true;
            }
        }

        bool has_mismatch = false;

        for(int p=1; p<=6; p++) {
            if(control_used[p]&&target_used[p]) {
                has_mismatch = true;
            }
        }

        if(is_in_c3 && has_mismatch) {
            bool flag = true;
            // if pi is in c3 but is not mismatch-free, we will test some options for a maximal abelian subgroup A of the pauli group
            // such that pi conjugates A into the pauli group. it will happen that at least one of these options always works
            if(flag) {
                // checks <x6, x5, x3x4, z3z4, z2, z1>
                if(yields_x(6,1,pi,pi_inv)&&yields_x(6,2,pi,pi_inv)&&yields_x(6,12,pi,pi_inv)&&yields_z(6,12,pi,pi_inv)&&yields_z(6,16,pi,pi_inv)&yields_z(6,32,pi,pi_inv)) {
                    flag = false;
                }
            }

            if(flag) {
                // checks <x6, x4x5, z4z5, z3, z2, z1>
                if(yields_x(6,1,pi,pi_inv)&&yields_x(6,6,pi,pi_inv)&&yields_z(6,6,pi,pi_inv)&&yields_z(6,8,pi,pi_inv)&&yields_z(6,16,pi,pi_inv)&yields_z(6,32,pi,pi_inv)) {
                    flag = false;
                }
            }

            if(flag) {
                // checks <x6, x3x5, z4, z3z5, z2, z1>
                if(yields_x(6,1,pi,pi_inv)&&yields_x(6,10,pi,pi_inv)&&yields_z(6,4,pi,pi_inv)&&yields_z(6,10,pi,pi_inv)&&yields_z(6,16,pi,pi_inv)&yields_z(6,32,pi,pi_inv)) {
                    flag = false;
                }
            }

            if(flag) {
                // checks <x6, x3x4x5, z4z5, z3z5, z2, z1>
                if(yields_x(6,1,pi,pi_inv)&&yields_x(6,14,pi,pi_inv)&&yields_z(6,6,pi,pi_inv)&&yields_z(6,10,pi,pi_inv)&&yields_z(6,16,pi,pi_inv)&yields_z(6,32,pi,pi_inv)) {
                    flag = false;
                }
            }

            if(flag) {
                // checks <x6, z5, x3x4, z3z4, z2, z1>
                if(yields_x(6,1,pi,pi_inv)&&yields_z(6,2,pi,pi_inv)&&yields_x(6,12,pi,pi_inv)&&yields_z(6,12,pi,pi_inv)&&yields_z(6,16,pi,pi_inv)&yields_z(6,32,pi,pi_inv)) {
                    flag = false;
                }
            }

            if(flag) {
                print_bin(20,i);
                // nothing ever gets printed here, so we get the desired result
                cout << '\n';
            }
        }
    }

    cout << "search complete" << '\n';

    return 0;
}
