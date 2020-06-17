// msslbu.java - merge sort single link list bottom up
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

public class x {
    
    // merge two already sorted lists
    static Node temp = new Node();
    static Node merge(Node list0, Node list1) {
        if(list0 == null)
            return list1;
        if(list1 == null)
            return list0;
        Node dest = temp;
        while(true){
            if(list0.data <= list1.data){
                dest.next = list0;
                dest = list0;
                list0 = list0.next;
                if(list0 == null){
                    dest.next = list1;
                    break;
                }
            } else {
                dest.next = list1;
                dest = list1;
                list1 = list1.next;
                if(list1 == null){
                    dest.next = list0;
                    break;
                }
            }
        }
        return temp.next;
    }

    // bottom up merge sort for single link list
    static Node sortbu(Node head) {
        final int NUMLIST = 32;
        Node[] alist = new Node[NUMLIST];
        Node node;
        Node next;
        int i;
        // if < 2 nodes, return
        if(head == null || head.next == null)
            return null;
        node = head;
        // merge node into array
        while(node != null){
            next = node.next;
            node.next = null;
            for(i = 0; (i < NUMLIST) && (alist[i] != null); i++){
                node = merge(alist[i], node);
                alist[i] = null;
            }
            if(i == NUMLIST)   // don't go past end of array
                i--;
            alist[i] = node;
            node = next;
        }
        // node == null
        // merge array into single list
        for(i = 0; i < NUMLIST; i++)
            node = merge(alist[i], node);
        return node;
    }

    // test sort
    static Node testsort(Node head, int n)
    {
        Node node = head;
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
        head = sortbu(head);
        end = System.currentTimeMillis();
        System.out.println("milliseconds " + (end-bgn));
        // verify sort
        i = 1;
        node = head;
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
        return head;
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