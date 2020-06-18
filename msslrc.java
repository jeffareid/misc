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

    static void merge(Node prev, Node f0, int n0,
                       Node f1, int n1, NodePair np)
    {
        Node first = f0.data <= f1.data ? f0 : f1;
        Node last = prev;
        int i = 0;
        int j = 0;
        while(true){
            if(f0.data <= f1.data){         // if f0 < f1
                last.next = f0;             //  move f0
                last = f0;
                f0 = f0.next;
                if(++i < n0)                //  if not end run 0
                    continue;               //   continue back to while
                last.next = f1;             //  else link run 1
                while(++j < n1)
                    f1 = f1.next;
                last = f1;
                f1 = f1.next;
                break;
            } else {                
                last.next = f1;             //  move f1
                last = f1;
                f1 = f1.next;
                if(++j < n1)                //  if not end run 1
                    continue;               //   continue back to while
                last.next = f0;             //  else link run 0
                while(++i < n0)
                    f0 = f0.next;
                last = f0;
//              f0 = f0.next;
                break;
            }
        }
        np.first  = first;
        np.last   = last;
        last.next = f1;
    }

    static void mergesort(Node prev, Node f0, int n0, NodePair np)
    {
        if(n0 <= 1){
            np.first = np.last = f0;
            return;
        }
        Node f1;
        int n1; 
        n1 = n0; n0 >>= 1; n1 -= n0;
        mergesort(prev, f0, n0, np);
        f0 = np.first;
        f1 = np.last.next;
        mergesort(np.last, f1, n1, np);
        f1 = np.first;
        merge(prev, f0, n0, f1, n1, np);
    }

    // test sort
    static Node testsort(Node head, int n)
    {
        Node node = head;
        NodePair np= new NodePair();
        Node prev = new Node();
        prev.next = head;
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
        mergesort(prev, head, n, np);
        end = System.currentTimeMillis();
        System.out.println("milliseconds " + (end-bgn));
        // verify sort
        i = 1;
        node = prev.next;
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
        return prev.next;
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