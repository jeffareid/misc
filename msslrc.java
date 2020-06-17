// msslrc.java - merge sort single link list recursive
package x;
import java.util.Random;

class Node{
    Node next;
    int data;
    Node(){
        data = 0;
    }
    Node(int d){
        data = d;
    }
}

class NodePair{
Node first;
Node last;
}

public class x {

    static void merge(Node before, Node F1, int N1,
                       Node F2, int N2, NodePair NP)
    {
        Node first, last, temp;
        int I, J;
        first = last = F1.data <= F2.data ? F1 : F2;
        for (I=J=0; I < N1 || J < N2; ) {
            if (I < N1 && (J>=N2 || F1.data <= F2.data))
                 { temp = F1; F1 = F1.next; I++; }
            else { temp = F2; F2 = F2.next; J++; }
            last.next = temp;
            last = temp;
        }
        before.next = first;
        last.next = F2;
        NP.first = first;
        NP.last = last;       
    }

    static void mergesort(Node before, Node F1, int N1, NodePair NP)
    {
        if (N1 <= 1)
            NP.first = NP.last = F1;
        else {
            Node F2;
            int N2; 
            N2 = N1; N1 >>= 1; N2 -= N1;
            mergesort(before, F1, N1, NP);
            F1 = NP.first;
            F2 = NP.last.next;
            mergesort(NP.last, F2, N2, NP);
            F2 = NP.first;
            merge(before, F1, N1, F2, N2, NP);
        }
    }

    // test sort
    static Node testsort(Node head, int n)
    {
        Node node = head;
        NodePair np= new NodePair();
        Node before = new Node();
        before.next = head;
        long bgn, end;
        int i;
        // fill list with random data
        Random r = new Random();
        for(i = 0; i < n; i++){
            node.data = r.nextInt();
            node = node.next;
        }
        // time sort
        bgn = System.currentTimeMillis();
        mergesort(before, head, n, np);
        end = System.currentTimeMillis();
        System.out.println("milliseconds " + (end-bgn));
        // verify sort
        i = 1;
        node = before.next;
        while(node.next != null){
            if(node.data > node.next.data)
                break;
            node = node.next;
            i++;
        }
        if(i == n)
            System.out.println("sort passed");
        else
            System.out.println("sort failed");
        return before.next;
    }

    // main
    public static void main(String[] args) {
        // create list
        final int COUNT = 16*1024*1024; 
        Node head, node;
        int i;
        head = new Node();
        node = head;
        for(i = 1; i < COUNT; i++){
            node.next = new Node();
            node = node.next;
        }
        head = testsort(head, COUNT);   // test sort with sequential nodes
               testsort(head, COUNT);   // test sort with scattered  nodes
    }
}